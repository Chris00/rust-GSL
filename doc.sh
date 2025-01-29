RUSTDOCFLAGS="-Z unstable-options --cfg docsrs -D warnings --generate-link-to-definition --html-in-header src/_docs/header.html" cargo doc --no-deps --all-features
