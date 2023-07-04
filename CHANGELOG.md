# Changelog
All notable changes to this project will be documented in this file.  
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.6.1] - 2023-07-04

### Added
- The `json` output format of `Disteval` now includes the values of the integrals in addition to the values of the sums.
- [zlib](http://zlib.net/) version 1.2.13, needed by FORM.
- Example `nodist_examples/ggh`. Demonstrates the computation of the 1- and 2-loop amplitudes for Higgs production in gluon fusion.
- Example `nodist_examples/triangle2L_wu`. Integral that requires a rescaling of the Feynman parameters.

### Changed
- FORM is now configured using `--with-zlib`, to avoid [FORM issue 95](https://github.com/vermaseren/form/issues/95).
- If `lattice_candidates` is even and non-zero, use `lattice_candidates+1` candidates.
- The default `decomposition_method` of `make_package` and `code_writer.make_package` to `geometric_no_primary`.
- The default `decomposition_method` of `loop_integral.loop_package` to `geometric`.
- [GiNaC](https://www.ginac.de/) updated to 1.8.6+ (commit 4bc2092 from Jun 21 2023).

### Fixed
- Critical bug in `disteval` introduced in 1.6 when using the median QMC lattice rules. Bug led to incorrect results for integrals with a severely underestimated error after recomputing with a larger lattice.
- Parsing of the coefficient expressions with `ginsh`. Previously, mixing the exponentiation operator `^` and the unary `+` and/or `-` operators would result in `ginsh` misparsing the coefficients.

## [1.6] - 2023-05-29

### Added
- New integrator `Disteval`.
- Integration based on median QMC rule implemented, can be enabled with option `lattice_candidates`.
- Numerator support for expansion by regions.
- `suggested_extra_regulator_exponent`  function, returns a list of suggested extra regulators sufficient to regularise a loop integral.
- `extra_regulator_constraints` function, returns a dict of inequalities which must be obeyed by the extra regulators in order to regularise a loop integral.
- `form_memory_use` argument for  `loop_regions`, tailors `form.set` to use approximately the requested amount of memory.
- `form_threads` argument for `loop_regions`, the number of threads TFORM will use.
- `extra_regulator_name`, `extra_regulator_exponent` for `loop_regions`, the name to be used for the extra regulator required to obtain well defined integrals and a list of exponents of the extra regulator.
- Documentation for `prefactor` in `IntegralLibrary` output.
- `ginsh` binary built for GiNaC, now used for coefficient parsing.
- Example `examples/pentabox_offshell`. A 2-loop penta-box integral.
- Example `examples/hexatriangle`. A massive 2-loop qq->ttH master integral
- Example `examples/region_tools`. Demonstrates the standalone usage of `suggested_extra_regulator_exponent`, `extra_regulator_constraints` and `find_regions`.
- Added [jupyter](https://jupyter.org/) examples in `examples/jupyter`.

### Changed
- Vastly improved computation of the Newton polytopes, speeding up some steps required for expansion by regions and geometric sector decomposition.
- Prefactors now expanded with `ginsh` rather than SymPy (while still using the SymPy syntax).
- Coefficient parsing now relies on `ginsh`, allows much more general coefficients than the previously required `num/den` rational function form.
- `sum_package` accepts much more general coefficients, which can be provided simply as strings with arbitrary arithmetic expressions (to be parsed by `ginsh`).
- `sum_package` accepts sum coefficients as dictionaries of the form `{'sum name' : terms}`, where `terms` is either a list of coefficient expressions (one per integral), or a dictionary of the form `{integral_index : coefficient}`, allowing for sparse coefficient matrices.
- The default sector decomposition method in `loop_package` changed from `iterative` to `geometric`.
- Use `form` instead of `tform` and `form_threads=1` by default, parallelisation is provided by the build system instead.
- Disabled `ContinuationLines` in FORM output.
- Various scripts `export_sector`, `formwrapper`, `write_contour_deformation`, `write_integrand` moved to `pySecDecContrib`.
- `git_id` changed to `__commit__` to be more consistent with naming of other metadata (e.g. `__authors__` and `__version__`).
- Print `test.log` for failed high level tests.
- Require recent version of `numpy>=1.23`.
- Require recent version of `sympy>=1.10.1` and `sympy<1.11` (due to a bug in the sympy series expansion).
- Python binary can be set with `PYTHON` variable.
- Replaced python testing framework `nose` with `pytest`.
- [Catch2](https://github.com/catchorg/Catch2) version 3.3.2 is now included in `high_level_tests` and used as our C++ testing framework.
- [Normaliz](https://www.normaliz.uni-osnabrueck.de/) version 3.9.2 is now included in `pySecDecContrib`; `normaliz_executable` and `normaliz` arguments to `make_package` and other functions are now optional (and should probably not be used).
- [GiNaC](https://www.ginac.de/) updated to 1.8.4.
- [CLN](https://www.ginac.de/CLN/) updated to 1.3.6-b4d44895.
- [FORM](https://github.com/vermaseren/form) updated to 4.3.0.
- [Cuba](https://feynarts.de/cuba/) updated to 4.2.2.
- [Nauty and Traces](https://pallini.di.uniroma1.it) updated to 2.8.6.
- [GSL](https://www.gnu.org/software/gsl/) updated to 2.7.1.

### Removed
- Support for Python versions 3.6 and 3.7.
- A C++17 compliant compiler is now required (previously C++14 was sufficient).
- `add_monomial_reglator_power` argument for `loop_regions`, replaced by `extra_regulator_name` and `extra_regulator_exponent`.

### Fixed
- GPU support for CUDA 11, removed incorrect use of `shared_ptr` in device functions.
- `geometric_ku` now correctly handles 0-dimensional cones.
- Handling of the imaginary unit `i_` when they appear e.g. in user-provided polynomials.
- Expansion by regions for cases where the resulting expansion is trivial.
- Provide more useful error messages in `polytope` class, relevant when using expansion by regions or geometric sector decomposition methods.
- Deprecation warnings emitted by SymPy due to calls of type `sympify(str)`
- The Cuba examples are not built during installation

## [1.5.6] - 2022-11-15

### Fixed
- Critical bug in `expand_singular` introduced in alpha_v0.1, which could lead to subtly incorrect analytic and numerical results. The series expansion of rational functions was incorrectly truncated if the expansion happened to be zero at some order in the expansion (Thanks to Christoph Greub).

## [1.5.5] - 2022-09-12

### Fixed
- Release version number

## [1.5.4] - 2022-09-10

### Fixed
- Critical bug in `_make_FORM_function_definition` introduced in v1.5, which could lead to subtly incorrect numerical results when polynomials with one term and a large coefficient (more than 10^6 characters) were encountered in intermediate stages of the decomposition (Thanks to Christoph Greub).

## [1.5.3] - 2022-01-26

### Added
- Example `nodist_examples/gggH1L`. Higgs plus jet production at 1-loop. (Thanks to @CpSquared)
- Example `nodist_examples/banana4L`. A 4-loop banana integral.

### Changed
- Vastly improved parallelisation of code generation and compilation.
- Improved performance of `remap_one_to_zero` (used when `split=True`).
- Default integrator in `integral_interface` is now set on call rather than on instantiation. This allows `integral_interface` to be used even when the default qmc transform is not available.
- Example `examples/two_regulators` now demonstrates the use of asymmetric Korobov transforms.

### Fixed
- Code generated by `sum_package` when passing a mix of real and complex integrals. (Thanks to @pguzowski)
- Missing CUDA decorators to complex power function. (Thanks to @pguzowski)
- Prefactor in quoted analytic result for  `examples/elliptic2L_euclidean`.
- Building of documentation.
- Sympy deprecation warning when calling `make_regions` directly.

## [1.5.2] - 2021-08-24

### Removed
- Example `nodist_examples/ExpansionByRegions`.

### Fixed
- Building of the [online documentation](https://secdec.readthedocs.io).

## [1.5.1] - 2021-08-24

### Changed
- Print more detailed information regarding the stage of the calculation (generating integrals/ optimising contour/ summing integrals) when `verbose=1`.

### Fixed
- Issue #7: Bug which modified exponents in the original integrand expression when computing derivatives. This led to incorrect results for some integrals with very simple `F` polynomials. (Thanks to @apik for reporting this)

## [1.5] - 2021-08-09 (Expansion by Regions / Amplitudes)

### Added
- `sum_package` function, generates a library for integrating a sum (or multiple sums) of integrals/sectors. Compared to the libraries made using the old `make_package`, the `sum_package`-derived libraries:
  - automatically determine to which precision to evaluate each integral in a sum to achieve the overall targe precision in the fastest way possible;
  - automatically adjust contour deformation parameters, so the users are not required to tweak `deformation_parameters_maximum` to avoid sign check errors.
- `make_regions` function, performs expansion by regions on generic regulated integrals.
- `loop_regions` function, performs expansion by regions on loop integrals.
- `Coefficient` class, represents a rational function that can be passed to `sum_package`.
- `series_to_ginac`, `series_to_sympy`, `series_to_mathematica`, and `series_to_maple` functions added to `pySecDec.integral_interface`, they convert the result returned by the python library to a format compatible with various CAS programs. See [issue #3](https://github.com/gudrunhe/secdec/issues/3).
- `pylink_qmc_transforms` argument for `sum_package`, `make_package`, and `loop_package`, allows a list of QMC integral transforms that should be generated for the python library to be specified (e.g. `pylink_qmc_transforms=['korobov2x3','korobov4','sidi4']`). Default is now just `korobov3`.
- `form_memory_use` argument for `make_package` and `loop_package`, the maximum memory use allowed by FORM.
- `form_threads` argument for `make_package` and `loop_package`, the number of allowed TFORM threads.
- `DEBUG=1` flag for C++ library, produces a library suitable for debuggers and applies AddressSanitizer.
- Examples `easy_sum` and `yyyy1L`: demonstrates `sum_package` with generic integrals.
- Example `yyyy1L`: demonstrates `sum_package` with loop integrals.
- Example `make_regions_ebr`: demonstrates `make_regions`.
- Examples `box1L_ebr`, `bubble1L_dotted_ebr`,  `bubble1L_ebr`, `bubble2L_largem_ebr`, `bubble2L_smallm_ebr`, `formfactor1L_ebr`, `make_regions_ebr`, and `triangle2L_ebr`: demonstrate `loop_regions`.
- "Getting Started" sections for the documentation of `sum_package`, `make_regions`, and `loop_regions`.
- Documentation for amplitude header.
- [GiNaC](https://www.ginac.de/) 1.8.0 and [CLN](https://www.ginac.de/CLN/) 1.3.6 are now included in the distribution; libraries produced with `sum_package` depend on both.
- All the user-facing functions and classes are now exported directly from the `pySecDec` module, and are accessible as `pySecDec.<name>`.

### Changed
- Installation is now performed via `python3 -m pip install --user --upgrade pySecDec`.
- The C++ frontend has been partly rewritten, see the generated `integrate_<name>.cpp` for the new usage.
- `make_package` and `loop_package` now default to using `sum_package` underneath, improving their performance and robustness.
  - The old behaviour of `make_package` can still be obtained by calling `pySecDec.code_writer.make_package`.
  - The old behaviour of `loop_package` can still be obtained by passing the argument `package_generator=pySecDec.code_writer.make_package` to the function.
- Multiple CUDA architectures can now be specified when building the C++ library by setting the environment variable `SECDEC_WITH_CUDA_FLAGS` to `-gencode arch=compute_XX, code=sm_XX -gencode arch=compute_YY,code=sm_YY`. Please refer to the [CUDA NVCC documentation](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/) for details of this syntax; the script `print-cuda-arch.sh` in `examples/easy` may be useful for determining your architecture.
- `IntegralLibrary` supports evaluating weighted sums of integrals, extra optional arguments for controlling the evaluation have been added (see the [online documentation](https://secdec.readthedocs.io/en/stable/full_reference.html#pySecDec.integral_interface.IntegralLibrary) for details).
- To save compilation time, not every QMC integral transform is compiled for the python interface by default. Pass the `pylink_qmc_transforms` argument to `*_package` when generating the library to enable additional transformations.
- The [formset](https://github.com/tueda/formset) script is now used to tune the size of the FORM buffers to allow maximum resources without exceeding the user specified memory bound; the `form_memory_use` parameter to `*_package` sets the desired maximal memory usage.
- The `WorkSpace` parameter of FORM is automatically increased during library compilation. Users are no longer required to manually adjust `form_work_space` to avoid FORM-related failures.
- `make_package` uses less RAM.
- Various optimisations when computing derivatives and Jacobians (should slightly improve generate performance).
- Parallelisation in Makefile for C++ library improved.
- Example `userdefined_cpp` made compatible with the new C++ backend.
- [Cuba](http://www.feynarts.de/cuba/) was updated to 4.2.1 (28 Jun 2021).
- [Qmc](https://github.com/mppmu/qmc) was updated to 1.0.6.
- Integration routines now print logs to the standard error instead of the standard output. To ignore the distinction and save both, invoke your integration script as `python3 integrate.py >output.txt 2>&1` from the shell.

### Removed
- Python versions below 3.6 are no longer supported.
- A C++14 compliant compiler is now required (previously C++11 was sufficient).
- Example `easy_cuda` was removed: CUDA is now used by default when possible, so it is equivalent to the `easy` example.
- Environment variable `SECDEC_CONTRIB` is no longer needed.

### Deprecated
- The `SECDEC_WITH_CUDA=sm_XX` environment variable has been deprecated, use `SECDEC_WITH_CUDA_FLAGS=-arch=sm_XX` instead.
- The `requested_order=x` argument of `loop_package` is deprecated in favor of `requested_orders=[x]`.
- The `regulator=x` argument of `LoopIntegralFromGraph` and `LoopIntegralFromPropagators` is deprecated in favor of `regulators=[x]`.

### Fixed
- Polynomial to power 0 (`poly**0`) was not simplified to 1.
- 1-loop tadpole with geometric decomposition.
- Missing `if name == "__main__"` guards in some examples, required by `multiprocessing`.

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

## [1.2] - 2017-08-09 [YANKED]
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

