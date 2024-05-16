# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

pV <- packageVersion("MBBEFDLite")
myOtherPkgs <- c("Delaporte", "lamW", "Pade", "revss", "minimaxApprox")


# Test CITATION has most recent package version
expect_true(any(grepl(pV, toBibtex(citation("MBBEFDLite")), fixed = TRUE)))

# Test NEWS has most recent package version
expect_true(any(grepl(pV, news(package = "MBBEFDLite"), fixed = TRUE)))

# Test that CITATION doesn't contain the name of any other of my packages
expect_false(any(sapply(myOtherPkgs, grepl,
                       x = toBibtex(citation("MBBEFDLite"), fixed = TRUE))))

# Test that NEWS doesn't contain the name of any other of my packages
expect_false(any(sapply(myOtherPkgs, grepl,
                        x = news(package = "MBBEFDLite"), fixed = TRUE)))
