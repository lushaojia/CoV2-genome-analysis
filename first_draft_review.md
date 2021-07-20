# Review of Analysis and Code repos

**date**: 2021-07-19

## Review Notes

- See instructions in Pull Request for incorporating changes
- Comments in julia code are preceded with `#`
- Comments in markdown look like this: `<!-- this is a comment -->`
- Comments that require action begin with `TODO`, others are just for information
- You may delete comments inside other files once you've addressed them
- Please keep this file in the repository - you may add your own responses if aplicable.

## Analysis plan

- Some [sembr](http://sembr.org)
- other formatting

## Code repo

https://github.com/lushaojia/BioinformaticsBISC195.jl

- Many functions need tests
- Current tests don't pass:

  ```
  Test Summary:          | Pass  Fail  Error  Total
  BioinformaticsBISC195  |   23     3     10     36
    Using Strings        |   23     3     10     36
      normalizeDNA       |    5                   5
      composition        |          1     10     11
      gc_content         |    2     2             4
      complement         |    4                   4
      reverse_complement |    4                   4
      parse_fasta        |    8                   8
  ERROR: LoadError: Some tests did not pass: 23 passed, 3 failed, 10 errored, 0 broken.
  ```

## Analysis repo

Could use slightly more narrative, and I think maybe just sentences
instead of all bullet points would be good but overall really nice!
Really well done on dotplot!

Looks like you only have one of the analyses done though.