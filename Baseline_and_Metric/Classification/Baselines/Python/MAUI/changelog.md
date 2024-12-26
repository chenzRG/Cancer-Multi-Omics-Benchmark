# Changelog

## [unreleased] - ????-??-??

### Fixed

- Fixed a bug in fine-tuning models that were pre-trained, when the pre-trained model is first deserialized from disk.

## [0.1.7] - 2019-10-08

- `maui_model.transform()` no longer automatically computes feature correlations. Added a method `maui_model.get_feature_correlations()` to use when this is desirable. This is a huge performance enhancement when feature correlations aren't necessary.
