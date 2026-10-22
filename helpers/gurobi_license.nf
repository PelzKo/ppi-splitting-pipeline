// Resolve params.gurobi_license to a file, or [] when there is none to stage.
//
// Functions-only, like mqc_labels.nf: the caller passes the param in, so nothing
// here reads params. Both SPLIT_POSITIVES and SAMPLE_NEGATIVES need the same rule.
//
// The default value is the bare name "gurobi.lic", which used to go through
// checkIfExists unconditionally -- so a fresh clone with no license file failed
// at channel construction even for a kahip-only run, and -profile test only worked
// where an untracked gurobi.lic happened to sit in the launch dir. Now:
//   - unset, blank, or the true that `--gurobi_license ''` arrives as        -> []
//   - the default name, missing                                              -> [] plus a warning
//     (gurobipy then uses its own lookup: GRB_LICENSE_FILE, ~/gurobi.lic, or the
//     size-limited license bundled with the pip wheel)
//   - any other path, missing                                                -> error, since the
//     user named it explicitly
def gurobiLicenseFile(value) {
    if (!value || value instanceof Boolean || value.toString() == "true") {
        return []
    }
    def lic = file(value.toString())
    if (lic.exists()) {
        return lic
    }
    if (value.toString() == "gurobi.lic") {
        log.warn "No gurobi.lic in the launch directory -- Gurobi (if used) falls back to its own license lookup. Pass --gurobi_license to stage one."
        return []
    }
    error("--gurobi_license '${value}' does not exist.")
}
