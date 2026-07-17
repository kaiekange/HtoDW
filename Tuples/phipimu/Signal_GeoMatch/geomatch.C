// association vs geometric truth-matching, phipimu, gen level, PF only.
//
// association (pf): the association map + the pdgId/TwoBodyDecay ancestry chain,
// identical to H_phipimu except that its dR<0.05 / |dpt|<0.3 cuts are applied
// here rather than at fill time. The map identifies the particle, so those cuts
// are pure post-hoc quality and association stays scannable offline.
//
// geometric (geo): for each gen daughter, every same-charge charged-hadron PF
// passing BOTH dR<0.05 and |dpt|<0.3, then the smallest dR among the survivors.
// Both cuts are applied during selection, in the analyzer -- they decide WHICH
// candidate wins, so they cannot be deferred: a soft track that happens to be
// nearest must not block a good candidate slightly further out. geometric is
// therefore fixed at (0.05, 0.3) and a different working point needs a rerun.
//
// Run:  root -l -q 'geomatch.C("path/to/geo.root","GEOphipimuDsW")'
//
// Reading the cross-tag: for the PF geometric picked, <dau>_geo_assoc_key
// is that PF's OWN association key into prunedGenParticles.
//   key == Gen_<dau>_idx -> confirmed correct
//   key == null (-777)   -> the geo pick has no association at all
//   otherwise            -> the geo pick is truthfully a different particle
// The "confirmed" bucket is structurally EMPTY inside (geo_match && !pf_match):
// if the geo pick's association pointed at our daughter, association would
// have found that same PF and pf_match would be true. It is printed anyway as a
// self-check -- a nonzero count there means something is wrong.

#include <cstdio>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

namespace {
const char *kDau[3] = {"Kp", "Km", "pi"};

Long64_t N(TTree *t, const char *cut) { return t->GetEntries(cut); }

// mean/RMS of a branch under a cut, for the residual comparison
void Resid(TTree *t, const char *var, const char *cut, Long64_t &n, Double_t &mean, Double_t &rms)
{
	TH1F h("h_resid", "", 400, -2, 2);
	t->Draw(Form("%s>>h_resid", var), cut, "goff");
	n = (Long64_t)h.GetEntries(); mean = h.GetMean(); rms = h.GetRMS();
}
} // namespace

void geomatch(const char *file = "geo.root", const char *dir = "GEOphipimuDsW")
{
	TFile *f = TFile::Open(file);
	if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }
	TTree *t = (TTree *)f->Get(Form("%s/Events", dir));
	if (!t) { printf("no %s/Events in %s\n", dir, file); return; }
	const Long64_t n = t->GetEntries();
	printf("\n=== %s  (%s)   entries = %lld ===\n", file, dir, n);

	// ---- 1. the closure the whole cross-tag rests on -------------------------
	// The Association is built on prunedGenParticles (see
	// PackedCandidateGenAssociationProducer.cc), the same collection the gen
	// filter indexes, so association's key must equal Gen_<dau>_idx.
	printf("\n-- key closure (must be 0 mismatches) --\n");
	for (int i = 0; i < 3; i++)
	{
		Long64_t m = N(t, Form("%s_pf_match", kDau[i]));
		Long64_t ok = N(t, Form("%s_pf_match && %s_pf_key==Gen_%s_idx", kDau[i], kDau[i], kDau[i]));
		printf("   %-2s  matched=%-6lld  key==Gen_idx=%-6lld  mismatch=%lld\n", kDau[i], m, ok, m - ok);
	}

	// ---- 2. efficiency, both methods at the same working point ------------------
	// geometric already IS (0.05, 0.3) by construction, so the comparison to make is
	// against association cut to the same point. association's looser rows are shown only to
	// expose how little its own cuts bite -- hit-based matching lands at dR ~
	// 0.01, so its dR<0.05 never fires and only |dpt|<0.3 removes anything.
	printf("\n-- efficiency (geo is fixed at dR<0.05, |dpt|<0.3 by construction) --\n");
	printf("   %-3s %-26s %8s %8s %10s\n", "dau", "assoc cut at", "assoc", "geo", "geo-assoc");
	for (int i = 0; i < 3; i++)
	{
		const char *L = kDau[i];
		const char *names[3] = {"(nothing)", "dR<0.05", "dR<0.05, |dpt|<0.3  <= same"};
		TString ca[3] = {Form("%s_pf_match", L),
		                 Form("%s_pf_match && %s_pf_dR<0.05", L, L),
		                 Form("%s_pf_match && %s_pf_dR<0.05 && abs(%s_pf_dpt)<0.3", L, L, L)};
		TString cg = Form("%s_geo_match", L);
		Double_t b = 100. * N(t, cg) / n;
		for (int k = 0; k < 3; k++)
		{
			Double_t a = 100. * N(t, ca[k]) / n;
			printf("   %-3s %-26s %7.2f%% %7.2f%% %+9.2fpp\n", L, names[k], a, b, b - a);
		}
		printf("\n");
	}

	// ---- 3. where geometry gains, what is it actually gaining? ---------------
	printf("-- decomposition of the geo-only gain (geo_match && !pf_match) --\n");
	for (int i = 0; i < 3; i++)
	{
		const char *L = kDau[i];
		Long64_t gain = N(t, Form("%s_geo_match && !%s_pf_match", L, L));
		Long64_t conf = N(t, Form("%s_geo_match && !%s_pf_match && %s_geo_assoc_key==Gen_%s_idx", L, L, L, L));
		Long64_t null = N(t, Form("%s_geo_match && !%s_pf_match && %s_geo_assoc_key==-777", L, L, L));
		printf("   %-2s gain=%-5lld  confirmed=%-4lld (must be 0)  assoc-null=%-5lld  other-particle=%lld\n",
		       L, gain, conf, null, gain - conf - null);
	}

	// ---- 4. is the assoc-null bucket real tracks? ----------------------------
	// The decisive test, and it needs no association to referee it: a genuinely
	// recovered track has dpt ~ 0 like any matched track. Soft junk picked up
	// because the real track was never reconstructed sits at dpt ~ -1.
	printf("\n-- residual of each bucket: real track => dpt ~ 0 --\n");
	for (int i = 0; i < 3; i++)
	{
		const char *L = kDau[i];
		Long64_t nn; Double_t m, r;
		Resid(t, Form("%s_pf_dpt", L), Form("%s_pf_match", L), nn, m, r);
		printf("   %-2s assoc-matched (reference) N=%-5lld  mean=%+.3f  RMS=%.3f\n", L, nn, m, r);
		Resid(t, Form("%s_geo_dpt", L), Form("%s_geo_match && !%s_pf_match && %s_geo_assoc_key==-777", L, L, L), nn, m, r);
		printf("   %-2s geo, assoc-null bucket    N=%-5lld  mean=%+.3f  RMS=%.3f\n", L, nn, m, r);
		Resid(t, Form("%s_geo_dpt", L), Form("%s_geo_match && !%s_pf_match && %s_geo_assoc_key!=-777", L, L, L), nn, m, r);
		printf("   %-2s geo, other-particle       N=%-5lld  mean=%+.3f  RMS=%.3f\n\n", L, nn, m, r);
	}

	// ---- 5. what switching would cost ----------------------------------------
	// The same two numbers read both ways. Agreement is what certifies the geometric
	// code: it is new, with no closure test of its own, but a buggy loop could
	// not agree with hit-based truth ~99% of the time. Disagreement is the purity
	// cost of switching: truth says the association pick is right, so every one of
	// these is a fake geometry would introduce.
	printf("-- geo vs hit-based truth on the commonly-matched set --\n");
	for (int i = 0; i < 3; i++)
	{
		const char *L = kDau[i];
		Long64_t both = N(t, Form("%s_pf_match && %s_geo_match", L, L));
		Long64_t same = N(t, Form("%s_pf_match && %s_geo_match && %s_pf_idx==%s_geo_idx", L, L, L, L));
		printf("   %-2s both=%-6lld  same PF=%-6lld (%.2f%%, validates geometric)  different=%-5lld (%.2f%%, purity cost)\n",
		       L, both, same, both ? 100. * same / both : 0., both - same, both ? 100. * (both - same) / both : 0.);
	}

	// K+ and pi carry the same charge, so geometry can hand ONE PF to both daughters.
	// The association cannot do this: it matches on truth identity.
	printf("\n-- same-sign K+/pi collision (a geo-only failure mode) --\n");
	printf("   Kp_geo_idx==pi_geo_idx : %lld / %lld  (%.2f%%)\n",
	       N(t, "Kp_geo_match && pi_geo_match && Kp_geo_idx==pi_geo_idx"), n,
	       100. * N(t, "Kp_geo_match && pi_geo_match && Kp_geo_idx==pi_geo_idx") / n);

	// ---- 6. the number the decision turns on ---------------------------------
	// Per-daughter gains compound as ~x^3, and so do the fakes, so judge on the triplet.
	printf("\n-- triplet (all three daughters) at the current working point --\n");
	// geometric carries its cuts already; re-applying them here would be a no-op that
	// invites the reader to think the two methods are cut differently.
	TString ca = "", cg = "";
	for (int i = 0; i < 3; i++)
	{
		const char *L = kDau[i];
		if (i) { ca += " && "; cg += " && "; }
		ca += Form("%s_pf_match && %s_pf_dR<0.05 && abs(%s_pf_dpt)<0.3", L, L, L);
		cg += Form("%s_geo_match", L);
	}
	TString cg_distinct = cg + " && Kp_geo_idx!=pi_geo_idx && Kp_geo_idx!=Km_geo_idx && Km_geo_idx!=pi_geo_idx";
	Long64_t a = N(t, ca), b = N(t, cg), c = N(t, cg_distinct);
	printf("   assoc                     : %-6lld / %lld = %.2f%%\n", a, n, 100. * a / n);
	printf("   geo                       : %-6lld / %lld = %.2f%%\n", b, n, 100. * b / n);
	printf("   geo, distinct tracks only : %-6lld / %lld = %.2f%%\n", c, n, 100. * c / n);

	// The arms run over the SAME events, so the two rates are correlated and the
	// naive difference has no meaningful binomial error. Only the discordant
	// triplets carry information (McNemar): quote the difference against
	// sqrt(discordant), or a fraction of a pp will look like a real effect.
	Long64_t only_a = N(t, Form("(%s) && !(%s)", ca.Data(), cg_distinct.Data()));
	Long64_t only_g = N(t, Form("!(%s) && (%s)", ca.Data(), cg_distinct.Data()));
	Double_t d = 100. * (Double_t)(only_g - only_a) / n;
	Double_t e = 100. * sqrt((Double_t)(only_a + only_g)) / n;
	Double_t z = (only_a + only_g) ? (only_g - only_a) / sqrt((Double_t)(only_a + only_g)) : 0.;
	printf("   discordant: assoc-only=%lld  geo-only=%lld\n", only_a, only_g);
	printf("   => geo - assoc = %+.2f +/- %.2f pp   (z = %+.2f)%s\n\n",
	       d, e, z, fabs(z) < 2 ? "   [consistent with zero: the two methods are equivalent]" : "");
}
