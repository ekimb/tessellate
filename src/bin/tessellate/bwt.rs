pub fn bwt(input: &str, output: &str, rate: usize) {
    tessellate::fmd::build::build_file(input, output, rate)
        .unwrap_or_else(|e| panic!("cannot build {output}.tbwt: {e}"));
}
