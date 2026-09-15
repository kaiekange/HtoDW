# JERC Application Tutorial

---

## Instruction to run the code is on the [JERC website](https://cms-jerc.web.cern.ch/ApplicationTutorial/#setup)


---

## To-Do List

* **Upgrade to NanoAOD v15**, ensuring compatibility with the custom `jetId` definition.
* **Revisit dR-based GenJet matching for Type-1 MET JER smearing**, since `CorrT1METJet` does not provide a `genJetIdx` branch.
* **Synchronize correctionlib versions between C++ and Python**, to guarantee identical `hashprng` behavior and eliminate discrepancies in stochastic JER smearing.
* [ApplyOnMiniAOD]: tweak JecAppliedJetProducer.cc for AK8 and illustrate the application for different systematics and years

---

## Code changes

* **[09-Sep-2026]** (Update/Add): Included latest tags for 2025/2026

* **[28-Aug-2026]** (Add): Included standalone corrlib 

* **[05-June-2026]** (Update): Synchronised to the [updated JER inputs](https://cms-talk.web.cern.ch/t/new-jer-smearing-inputs-available-for-2024-and-2025/145723) 

* **[24-April-2026]** (Update): Synchronised to the updated run-dependent residual JECs for all years of Run-2

* **[15-April-2026]** (Update): Synchronised to the updated run-dependent residual JECs for all years of Run-3 
 

* **[26-Jan-2026]** (BUGFIX): The JES nominal correction in C++ is now ptAfter/ptStart instead of ptAfter/ptRaw as it is applied to NanoAOD.Jet_pt

