import numpy as np
import awkward as ak
import correctionlib
import correctionlib.schemav2 as cs
import hist

# Jet flavor labels used in correctionlib efficiency maps
_FLAVORS = [0, 4, 5]  # udsg, c, b
_FLAVOR_LABEL = {0: "udsg", 4: "c", 5: "b"}

# Binning used for efficiency histograms
_ABSETA_EDGES = [0, 1.305, 2.5]
_PT_EDGES = [20., 30., 50., 70., 100., 140., 200., 300., 600., 1000.]

# Dictionary mapping NanoAOD btag algorithm names to the corresponding keys in the correctionlib CorrectionSet
correctionlib_keys = {
    "btagUParTAK4B" : "UParTAK4",
    "btagPNetB" : "particleNet",
    "btagRobustParTAK4B" : "robustParticleTransformer",
    "btagDeepFlavB" : "deepJet"
}

def compute_btagging_efficiency(
        jets: ak.Array,
        btag_algo: str | list,
        working_points: dict,
    ) -> correctionlib.CorrectionSet:
    '''
    Core function: computes the per-jet btagging efficiency for a given set of working points.
    For each WP, the efficiency is computed as a function the jet flavor, absolute eta, and pt.
    The jets array is used to fill a 2D histogram of (abseta, pt) for each flavor category (b, c, light).
    The efficiency is then extracted as the ratio of jets passing the WP to the total number of jets in each (abseta, pt) bin.
    The function returns a correctionlib.CorrectionSet object that can be used to evaluate the efficiency for any given jet flavor, abseta, and pt.
    The CorrectionSet will contain one correction key for each tagger in btag_algo (if it's a list) or for the single tagger (if it's a string).
    The axes of the correction will be "wp", "flav", "abseta", and "pt".
    The signature of the correction should be efficiency_map.evaluate(wp, jetflav, jetabseta, jetpt)

    Args:
        jets: awkward array of jet collection, must contain fields "pt", "eta", "hadronFlavour", and the btag discriminator(s) specified in btag_algo, as well
            as boolean pass mask fields for each WP if working_points is a list (e.g. "btagUParTAK4B_L" for the "L" WP of the "btagUParTAK4B" discriminator).
        btag_algo: str or list of str, NanoAOD branch name(s) for the btag discriminator(s), e.g. "btagUParTAK4B" or ["btagUParTAK4B", "btagPNetB"]
        working_points: dict mapping of WP name (e.g. "L", "M", "T") to discriminator threshold value

    Returns:
        correctionlib.CorrectionSet object containing the efficiency maps for the specified btag algorithm(s) and working points.
    '''
    if isinstance(btag_algo, str):
        btag_algos = [btag_algo]
    else:
        btag_algos = list(btag_algo)

    # Sort working_points by increasing threshold (loosest to tightest)
    working_points = dict(sorted(working_points.items(), key=lambda item: item[1]))

    corrections = []

    for algo in btag_algos:
        # Flatten jet quantities
        jetpt_flat = ak.to_numpy(ak.flatten(jets["pt"]))
        jetabseta_flat = np.abs(ak.to_numpy(ak.flatten(jets["eta"])))
        jetflav_flat = ak.to_numpy(ak.flatten(jets["hadronFlavour"]))
        jetdisc_flat = ak.to_numpy(ak.flatten(jets[algo]))

        # For each WP build efficiency as a function of flavor, abseta, pt
        # Structure: Category[wp] -> Category[flav] -> Binning[abseta] -> Binning[pt] -> float
        wp_content = []
        for wp in working_points:
            # Practical design choice (consistent with the docstring):
            # The `working_points` parameter is a dict of WP (e.g. "L", "M", "T") to discriminator threshold value.

            if isinstance(working_points, dict):
                threshold = working_points[wp]
                pass_mask = jetdisc_flat > threshold
            else:
                raise ValueError(
                    f"Cannot determine WP threshold for '{wp}': working_points must be a dict "
                    f"mapping WP name -> threshold, or jets must contain a field '{algo}_{wp}'."
                )

            flav_content = []
            for flav in _FLAVORS:
                flav_mask = jetflav_flat == flav

                jetpt_f = jetpt_flat[flav_mask]
                jetabseta_f = jetabseta_flat[flav_mask]
                pass_f = pass_mask[flav_mask]

                # Build 2D histograms: total and passing
                h_total = hist.Hist(
                    hist.axis.Variable(_ABSETA_EDGES, name="abseta"),
                    hist.axis.Variable(_PT_EDGES, name="pt"),
                )
                h_pass = hist.Hist(
                    hist.axis.Variable(_ABSETA_EDGES, name="abseta"),
                    hist.axis.Variable(_PT_EDGES, name="pt"),
                )
                h_total.fill(abseta=jetabseta_f, pt=jetpt_f)
                h_pass.fill(abseta=jetabseta_f[pass_f], pt=jetpt_f[pass_f])

                total_vals = h_total.values()  # shape (n_eta, n_pt)
                pass_vals = h_pass.values()

                # Efficiency: avoid division by zero; bins with no jets get efficiency 0
                with np.errstate(divide="ignore", invalid="ignore"):
                    eff_vals = np.where(total_vals > 0, pass_vals / total_vals, 0.0)

                # Build innermost pt Binning nodes for each abseta bin
                abseta_content = []
                for ieta in range(len(_ABSETA_EDGES) - 1):
                    pt_bin_content = [float(eff_vals[ieta, ipt]) for ipt in range(len(_PT_EDGES) - 1)]
                    pt_binning = cs.Binning.parse_obj({
                        "nodetype": "binning",
                        "input": "pt",
                        "edges": list(_PT_EDGES),
                        "content": pt_bin_content,
                        "flow": "clamp",
                    })
                    abseta_content.append(pt_binning)

                eta_binning = cs.Binning.parse_obj({
                    "nodetype": "binning",
                    "input": "abseta",
                    "edges": list(_ABSETA_EDGES),
                    "content": abseta_content,
                    "flow": "clamp",
                })

                flav_content.append({"key": int(flav), "value": eta_binning})

            wp_content.append({
                "key": wp,
                "value": cs.Category.parse_obj({
                    "nodetype": "category",
                    "input": "flav",
                    "content": flav_content,
                })
            })

        correction = cs.Correction.parse_obj({
            "version": 1,
            "name": f"{correctionlib_keys.get(algo, algo)}_efficiency",
            "description": (
                f"MC b-tagging efficiency for {algo}. "
                "Axes: wp (string), flav (int: 0=udsg, 4=c, 5=b), "
                "abseta (float), pt (float)."
            ),
            "inputs": [
                {"name": "wp",      "type": "string",  "description": "Working point: L, M, T, XT, XXT"},
                {"name": "flav",    "type": "int",     "description": "Jet hadron flavor: 0=udsg, 4=c, 5=b"},
                {"name": "abseta",  "type": "real",    "description": "Jet |eta|"},
                {"name": "pt",      "type": "real",    "description": "Jet pt [GeV]"},
            ],
            "output": {"name": "efficiency", "type": "real", "description": "MC b-tagging efficiency"},
            "data": cs.Category.parse_obj({
                "nodetype": "category",
                "input": "wp",
                "content": wp_content,
            }),
        })
        corrections.append(correction)

    cset = cs.CorrectionSet.parse_obj({
        "schema_version": 2,
        "description": "MC b-tagging efficiency maps",
        "corrections": corrections,
    })
    # Return the compiled (evaluable) CorrectionSet
    return correctionlib.CorrectionSet.from_string(cset.model_dump_json(exclude_unset=True))

def validate_efficiency_map(efficiency_map : correctionlib.Correction) -> None:
    '''
    Validates that efficiency_map is compatible with the expected call signature:
        efficiency_map.evaluate(wp, jetflav, jetabseta, jetpt)
    where wp is a string, jetflav is an int, and jetabseta and jetpt are floats (reals).

    Args:
        efficiency_map: object to validate, expected to be a correctionlib.Correction
            with inputs [string, int, real, real].

    Raises:
        TypeError: if efficiency_map has no callable evaluate() method or no inputs attribute.
        ValueError: if the number or types of inputs are incompatible with the expected signature.
    '''
    if not callable(getattr(efficiency_map, "evaluate", None)):
        raise TypeError(
            "efficiency_map must have a callable 'evaluate' method. "
            f"Got object of type {type(efficiency_map).__name__!r} with no such method."
        )

    inputs = getattr(efficiency_map, "inputs", None)
    if inputs is None:
        raise TypeError(
            "efficiency_map has no 'inputs' attribute. "
            "Expected a correctionlib.Correction object."
        )

    expected_types = ["string", "int", "real", "real"]
    expected_names = ["wp", "flav", "abseta", "pt"]

    if len(inputs) != len(expected_types):
        name = getattr(efficiency_map, "name", "<unknown>")
        raise ValueError(
            f"efficiency_map '{name}' has {len(inputs)} input(s), "
            f"expected {len(expected_types)}. "
            f"Required signature: evaluate(wp: string, flav: int, abseta: real, pt: real)."
        )

    for i, (inp, exp_type, exp_name) in enumerate(zip(inputs, expected_types, expected_names)):
        if inp.type != exp_type:
            name = getattr(efficiency_map, "name", "<unknown>")
            raise ValueError(
                f"efficiency_map '{name}': input #{i} ('{inp.name}') "
                f"has type '{inp.type}', expected '{exp_type}' "
                f"(for positional argument '{exp_name}')."
            )

def compute_btag_event_weight_multi_wp(
        jets : ak.Array,
        btag_algo : str,
        sf_corr_set : correctionlib.CorrectionSet,
        efficiency_map : correctionlib.Correction,
        working_points : list,
        variation_heavy : str,
        variation_light : str,
        ignore_warnings : bool = False
    ) -> ak.Array:
    '''
    Core function: computes the per-event btag SF weight for an arbitrary number of working points.
    See BTV POG twiki: https://btv-wiki.docs.cern.ch/PerformanceCalibration/fixedWPSFRecommendations/#scale-factor-recommendations-for-event-reweighting

    The WPs must be ordered from loosest to tightest, e.g. ["L", "M", "T"] or ["L", "M", "T", "XT", "XXT"].
    The per-event weight is the product of per-jet factors, where each jet falls into one of N+1 bins:
      - Fails loosest WP:           factor = (1 - eps_0)
      - Passes WP[k] but fails WP[k+1] (intermediate bins k=0..N-2): factor = (eps_k - eps_{k+1})
      - Passes tightest WP:         factor = eps_{N-1}
    The SF weight is the ratio of the DATA product to the MC product.

    Args:
        jets: awkward array of jet collection
        btag_algo: str, NanoAOD branch name for the btag discriminator, e.g. "btagUParTAK4B"
        sf_corr_set: correctionlib.CorrectionSet containing the SF corrections.
            The heavy-flavor SF is extracted with key f"{btag_algo}_comb",
            the light-flavor SF with key f"{btag_algo}_light",
            and the WP threshold values with key f"{btag_algo}_wp_values".
        efficiency_map: correctionlib.Correction for MC btag efficiency.
            Call efficiency_map.evaluate(wp, jetflav, jetabseta, jetpt) directly.
        working_points: list of WP names ordered loosest to tightest, e.g. ["L", "M", "T"]
        variation_heavy: str, name of the heavy-flavor variation
        variation_light: str, name of the light-flavor variation
        ignore_warnings: bool, if True, ignore warnings about zero or negative efficiencies and SFs (for testing purposes)

    Returns:
        per-event btag SF weight array
    '''

    # Extract SF corrections and WP values from the CorrectionSet
    btag_algo_key = correctionlib_keys.get(btag_algo, None)
    if btag_algo_key is None:
        raise ValueError(f"Btag algorithm {btag_algo} not recognized! Available algorithms: {list(correctionlib_keys.keys())}")
    sf_heavy = sf_corr_set[f"{btag_algo_key}_comb"]
    sf_light = sf_corr_set[f"{btag_algo_key}_light"]
    wp_values_corr = sf_corr_set[f"{btag_algo_key}_wp_values"]
    wp_values = {wp: wp_values_corr.evaluate(wp) for wp in working_points}
    # Sort wp_values and working_points by increasing WP threshold (loosest to tightest)
    wp_values = dict(sorted(wp_values.items(), key=lambda item: item[1]))
    working_points = list(wp_values.keys())

    # Flatten jet kinematics
    jetpt = ak.flatten(jets["pt"])
    jetabseta = abs(ak.flatten(jets["eta"]))
    jetflav = ak.flatten(jets["hadronFlavour"])
    jetcounts = ak.num(jets["pt"])

    # Separate heavy and light flavors
    heavy_mask = jets["hadronFlavour"] > 3
    light_mask = jets["hadronFlavour"] < 4

    jetpt_heavy = ak.flatten(jets["pt"][heavy_mask])
    jetabseta_heavy = abs(ak.flatten(jets["eta"][heavy_mask]))
    jetflav_heavy = ak.flatten(jets["hadronFlavour"][heavy_mask])

    jetpt_light = ak.flatten(jets["pt"][light_mask])
    jetabseta_light = abs(ak.flatten(jets["eta"][light_mask]))
    jetflav_light = ak.flatten(jets["hadronFlavour"][light_mask])

    # Per-WP: compute MC efficiency, SF, DATA efficiency
    effi_MC = {}
    effi_DATA = {}

    # Validate that the efficiency_map has the expected signature
    validate_efficiency_map(efficiency_map)

    for i, wp in enumerate(working_points):
        effi_MC[wp] = ak.unflatten(
            efficiency_map.evaluate(wp, jetflav, jetabseta, jetpt),
            counts=jetcounts
        )

        sf_light_flat = sf_light.evaluate(
            variation_light, wp, jetflav_light, jetabseta_light, jetpt_light
        )
        sf_heavy_flat = sf_heavy.evaluate(
            variation_heavy, wp, jetflav_heavy, jetabseta_heavy, jetpt_heavy
        )
        sf_flat = ak.to_numpy(ak.copy(jetpt))
        sf_flat[jetflav > 3] = sf_heavy_flat
        sf_flat[jetflav < 4] = sf_light_flat
        sf_DATA_MC = ak.unflatten(sf_flat, jetcounts)

        effi_DATA[wp] = sf_DATA_MC * effi_MC[wp]

    # Build per-jet factors for MC and DATA, then compute per-event product.
    # For N WPs ordered loosest (working_points[0]) to tightest (working_points[-1]):
    #   - Jets failing loosest WP:             factor = (1 - eps_0)
    #   - Jets in intermediate bin k (k=0..N-2),
    #     i.e. passing wp[k] but failing wp[k+1]: factor = (eps_k - eps_{k+1})
    #   - Jets passing tightest WP:            factor = eps_{N-1}

    disc = jets[btag_algo]
    pass_wp = {wp: disc > wp_values[wp] for wp in working_points}

    # Collect per-jet factors as list of jagged arrays, then concatenate
    p_MC_terms = []
    p_DATA_terms = []

    # 1) Jets passing the tightest WP
    tightest = working_points[-1]
    mask_tightest = pass_wp[tightest]
    p_MC_terms.append(effi_MC[tightest][mask_tightest])
    p_DATA_terms.append(effi_DATA[tightest][mask_tightest])

    # 2) Jets failing the loosest WP
    loosest = working_points[0]
    mask_fail_loosest = ~pass_wp[loosest]
    p_MC_terms.append((1 - effi_MC[loosest])[mask_fail_loosest])
    p_DATA_terms.append((1 - effi_DATA[loosest])[mask_fail_loosest])

    # 3) Intermediate bins: jets passing wp[k] but failing wp[k+1]
    for k in range(len(working_points) - 1):
        wp_lo = working_points[k]
        wp_hi = working_points[k + 1]
        in_bin = pass_wp[wp_lo] & ~pass_wp[wp_hi]
        p_MC_terms.append((effi_MC[wp_lo] - effi_MC[wp_hi])[in_bin])
        p_DATA_terms.append((effi_DATA[wp_lo] - effi_DATA[wp_hi])[in_bin])

    # 4) Concatenate all jet factors and extract numerator and denominator
    p_MC = ak.concatenate(p_MC_terms, axis=1)
    p_DATA = ak.concatenate(p_DATA_terms, axis=1)
    p_MC = ak.prod(p_MC, axis=-1)
    p_DATA = ak.prod(p_DATA, axis=-1)

    # 5) Compute SF as ratio of DATA product to MC product
    btag_sf = np.divide(p_DATA, p_MC)

    # Check if the numerator or the denominator are negative
    if not ignore_warnings:
        assert(ak.all(p_MC >= 0)), "Negative MC efficiency product encountered in multiple WP computation! Please check the SFs and efficiencies for the given variation. Heavy variation: {}, Light variation: {}".format(variation_heavy, variation_light)
        assert(ak.all(p_DATA >= 0)), "Negative DATA efficiency product encountered in multiple WP computation! Please check the SFs and efficiencies for the given variation. Heavy variation: {}, Light variation: {}".format(variation_heavy, variation_light)

    ##################################################################################################
    # Here we can add custom code to handle the case of vanishing or negative event weights, e.g. by replacing them with 1 (no correction).
    # For now we just raise an exception if that happens, to catch potential issues with the SFs and efficiencies.
    ##################################################################################################

    # Replace inf or NaN with 1
    finite_mask = np.isfinite(ak.to_numpy(btag_sf))
    btag_sf = ak.where(finite_mask, btag_sf, 1)

    # Check if the event weight is negative
    if not ignore_warnings:
        assert(ak.all(btag_sf >= 0)), "Negative btag SF encountered in multiple WP computation! Please check the SFs and efficiencies for the given variation. Heavy variation: {}, Light variation: {}".format(variation_heavy, variation_light)

    return btag_sf
