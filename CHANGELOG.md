# Changelog
All notable changes to this project will be documented in this file.  
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.5] - 2021-08-09 (Expansion by Regions / Amplitudes)

### Added
- `sum_package` function, generates a library for integrating a sum of integrals/sectors.
- `make_regions` function, performs expansion by regions on generic regulated integrals.
- `loop_regions` function, performs expansion by regions on loop integrals.
- `Coefficient` class, represents a rational function that can be passed to `sum_package`.
- `series_to_ginac`, `series_to_sympy`, `series_to_mathematica` and `series_to_maple` functions added to `pySecDec.integral_interface`, they convert the result returned by the python library to a format compatible with various CAS programs.
- `pylink_qmc_transforms` argument for `sum_package`, `make_package` and `loop_package`, allows a list of QMC integral transforms that should be generated for the python library to be specified (e.g. `pylink_qmc_transforms=['korobov2x3','korobov4','sidi4']`). Default is now just `korobov3`.
- `form_memory_use` and `form_threads` arguments for `make_package` and `loop_package`, the maximum memory use allowed by FORM and the number of TFORM threads, respectively. 
- `DEBUG=1` flag for c++ library, produces a library suitable for debuggers and applies AddressSanitizer.
- [example] `easy_sum` and `yyyy1L` (demonstrates `sum_package` with generic integrals).
- [example] `yyyy1L` (demonstrates `sum_package` with loop integrals).
- [example] `make_regions_ebr` (demonstrates `make_regions`).
- [example] `box1L_ebr`, `bubble1L_dotted_ebr`,  `bubble1L_ebr`, `bubble2L_largem_ebr`, `bubble2L_smallm_ebr`, `formfactor1L_ebr`, `make_regions_ebr`, `triangle2L_ebr` (demonstrate `loop_regions`).
- [docs] "Getting Started" sections for `sum_package`, `make_regions` and `loop_regions`.
- [docs] Documentation for amplitude header.

### Changed
- Installation is now performed via `python3 -m pip install --user --upgrade pySecDec`
- The c++ frontend has been partly rewritten, see the generated `integrate_<name>.cpp` for the new usage.
- Multiple CUDA architectures can now be specified when building the c++ library by setting the environment variable `SECDEC_WITH_CUDA_FLAGS="-gencode arch=compute_XX,
code=sm_XX -gencode arch=compute_YY,code=sm_YY"` (note: the script `print-cuda-arch.sh` in `examples/easy` may be useful for determining your architecture).
- `IntegralLibrary` supports evaluating weighted sums of integrals, extra optional arguments for controlling the evaluation have been added (see the online documentation for details).
- The [formset](https://github.com/tueda/formset) script is now used to tune the size of the FORM buffers to allow maximum resources without exceeding the user specified memory bound.
- `make_package` uses less RAM.
- Various optimisations when computing derivatives and Jacobians (should slightly improve generate performance).
- Parallelisation in Makefile for c++ library improved.
- Example `userdefined_cpp` is now compatible with the new c++ backend.
- [dist] update to Cuba-4.2.1 (28 Jun 2021).
- [dist] update to qmc-1.0.6

### Deprecated
- Python versions below 3.6 are no longer supported.
- A c++14 compliant compiler is now required (previously c++11 was sufficient)
- The behaviour of `make_package` and the generated c++ backend has changed. The old behaviour can be obtained by calling `pySecDec.code_writer.make_package`.
- The behaviour of `loop_package` and the generated c++ backend has changed. The old behaviour can be obtained by passing the argument `package_generator=pySecDec.code_writer.make_package` to the function.
- The `SECDEC_WITH_CUDA=sm_XX` environment variable has been deprecated, use `SECDEC_WITH_CUDA_FLAGS`.
- Not all QMC integral transforms are compiled for the python interface by default. Set `pylink_qmc_transforms` when generating the library.
- `requested_order=x` for `loop_package`, use `requested_orders=[x]`.
- `regulator=x` for `LoopIntegralFromGraph` and `LoopIntegralFromPropagators`, use `regulators=[x]`.
 
### Fixed
- Polynomial to power 0 (poly**0) was not simplified to 1.
- 1-loop tadpole with geometric decomposition.   
- Missing `if name == "__main__"` guards in some examples, required by `multiprocessing`.

### Removed
- Example `easy_cuda` (CUDA is now used by default where possible)


## [1.4.5] - 2020-12-13
- [examples] add non-planar 6-propagator box
- [dist] update to Cuba-4.2 (2020)
- [dist] update to qmc-1.0.4
- [make_package] init arrays to nan, avoids returning wrong results on GPU in rare cases
- [misc] fix sympify calls (avoids sympy deprecation warning)
- [high_level_tests] fix broken test selective_ibp
- [examples] fix compatibility with python 3.8.6
- [pylink] catch and print errors in python interface
- [algebra] improve use of .simplify()
- [make_package] optimize jacobian calculation
- [make_package] fix MaxDegreeFunction


## [1.4.4] - 2020-02-05
- [Vegas] increase default parameters `nstart` and `nincrease`
- [dist] update to qmc-1.0.3


## [1.4.3] - 2019-08-22
- [symmetries] add options to not run any symmetry finder
- [dist] update to form-4.2.1, qmc-1.0.2
- [doc] add faq section


## [1.4.2] - 2019-05-16
- [algebra] fix incomplete simplification of some expressions
- [doc] add instructions for sign_check_error
- [examples] correct HZ2L_nonplanar kinematics
- [tests] fix "regulator_in_powerlist" in combination with sympy-1.4


## [1.4.1] - 2018-11-29 
- [dist] update to qmc-1.0.1 (fix in PolySingular fit function)


## [1.4] - 2018-11-28
- [integral_interface] add the quasi-Monte Carlo (QMC) integrator which can optionally run on Graphics Processing Units (GPUs)
- [algebra] fix for sympy-1.3
- [dist] update to gsl-2.5

## [1.3.2] - 2018-08-02
- [prefactor expansion] fix bug if ``x``-expansion starts lower than ``1/x``
- [prefactor expansion] fix error if poles have multiple terms
- [subtraction] implement individual `ibp_power_goals` for the `indices` as suggested in [issue #2](https://github.com/gudrunhe/secdec/issues/2)
- [symmetry_finder] fix rare bugs that may occur on hash collisions
- [dreadnaut symmetry_finder] fix finding fake symmetries by illegal swappings; note that dreadnaut is disabled by default since version 1.2.2
- [dist] building the tarball with make's `-j` option is now supported
- [make_package] fix dropping of nonzero terms when using contour deformation in the presence of linear or higher poles

## [1.3.1] - 2018-04-24
- [integral_interface] support `MultiIntegrator` in python interface
- [make_package] compute determinants in parallel
- [make_package] fix "illegal instruction: 4" due to missing virtual destructors
- [dist] ship a newer version of [FORM](github.com/vermaseren/form/commit/77ee4eab218ff75bbc2f8e52a2d53efd06159fdf) with the [optimization bug](github.com/vermaseren/form/issues/272) fixed

## [1.3] - 2018-01-30
- [make_package] bugfixes concerning integrals with numerator in combination with higher than logarithmic poles
- [make_package] speed up algebraic part
- [util/integrator] implement "zero_border"
- [doc] list external dependencies and papers to cite

## [1.2.2] - 2017-08-25
- [loop_integral] fix issues with sympy-1.1.1
- [symmetry_finder] fix Pak's sorting algorithm
- [loop_package] fix error with regulator in 'powerlist'

## [1.2.1] - 2017-08-18
- fix release 1.2

## [1.2] - 2017-08-09
### Version 1.2 has a bug that can lead to wrong results without warning. Please use a different version.
- [make_package] more efficient algebra
- [util/integrator] add dedicated 1D integrator 'cquad'
- [util/integrator] implement 'MultiIntegrator' to choose an integrator depending on the dimension of the integrand

## [1.1.2] - 2017-06-26
- [make_package] fix unittest failing with python 3.6.1
- [util/integrator] fix one dimensional integration
- [util/integrand_container] fix memory access error in "complex_to_real"

## [1.1.1] - 2017-05-30
- [make_package]: drop lower bound on 'requested_order'
- [loop_integral]: fixed parameters of loop integral measure for integrals with both doubled and inverse propagators
- fix geometric_ku decomposition method

## [1.1] - 2017-05-20
- update documentation
- added example 'easy'
- update error propagation

## [1.0] - 2017-03-29 (Initial Release)
