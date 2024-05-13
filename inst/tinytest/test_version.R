# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

pV <- packageVersion("MBBEFDLite")

# Test CITATION has most recent package version
expect_true(any(grepl(pV, toBibtex(citation("MBBEFDLite")), fixed = TRUE)))

# Test NEWS has most recent package version
expect_true(any(grepl(pV, news(package = "MBBEFDLite"), fixed = TRUE)))
