
all: doc

doc:
	cargo doc --no-deps --all-features

deny:
	cargo deny check

audit:
	cargo audit

check:
	python3 checker.py

.PHONY: audit check deny doc
