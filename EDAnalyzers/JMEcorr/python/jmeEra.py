import FWCore.ParameterSet.Config as cms

# Single entry point for the JME era.
#
# "Year" is a knob on THREE separate modules -- the two JEC producers and the jet veto
# map filter -- and setting only some of them is silently wrong rather than an error:
# a mismatched veto map is just a different (eta,phi) map, with no run-number dimension
# to catch the mistake. (The JEC half is partly self-protecting on data, where
# L2L3Residual is binned in run number and throws on a mismatch, but not on MC.)
#
# Calling this function is the only supported way to set the era: it is impossible to
# update two modules out of three. The cfi files carry no usable default -- Year is ""
# there and each plugin throws on an empty value -- so a pset that forgets to call this
# fails at construction instead of quietly producing preEE numbers.
#
# Valid years are the top-level keys of JecConfigAK4.json / JvmConfig.json:
#   2016Pre 2016Post 2017 2018 2022Pre 2022Post 2023Pre 2023Post 2024 2025
# For data, era is the single data-era key under that year in JecConfigAK4.json
# ("Era2022PreAll" under 2022Pre, "Era2022PostAll" under 2022Post).


def applyJmeEra(process, year, isData=False, era=""):
    """Set the JME era consistently on jecAppliedJetProducer, jecAppliedMetProducer
    and jvmAppliedEventFilter.

    process : the cms.Process, with all three modules already loaded
    year    : JEC payload set + jet veto map, e.g. "2022Pre" / "2022Post"
    isData  : True enables the data JEC path (L2L3Residual, no JER smearing)
    era     : data-era key, required when isData is True, ignored otherwise
    """
    if not year:
        raise ValueError("applyJmeEra: 'year' is required")
    if isData and not era:
        raise ValueError("applyJmeEra: data (isData=True) also needs 'era', "
                         "e.g. era='Era2022PostAll' for year='2022Post'")
    if era and not isData:
        raise ValueError("applyJmeEra: 'era' is a data-only key; got era=%r with "
                         "isData=False" % era)

    for mod in (process.jecAppliedJetProducer, process.jecAppliedMetProducer):
        mod.Jets.Year = cms.string(year)
        mod.Jets.IsData = cms.bool(isData)
        mod.Jets.Era = cms.string(era)

    process.jvmAppliedEventFilter.Jets.Year = cms.string(year)
