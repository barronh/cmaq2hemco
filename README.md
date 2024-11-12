# cmaq2hemco

Utilities for converting CMAQ-ready emissions to HEMCO supported files.

The best way to learn about this module is the aws_mp2022.py script in the
examples directory. It relies on helper functions to pull data from AWS s3
buckets for the 2022 MP.

# To-Do

1. Currently stored in CB6 species. Needs HEMCO translation or added translation before save to disk.
2. Need to check that HEMCO can read these files (should be fine)