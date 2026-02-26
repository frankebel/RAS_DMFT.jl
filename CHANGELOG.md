# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- moments of complex function ([#164](https://github.com/frankebel/RAS_DMFT.jl/pull/164)) ([da721d2](https://github.com/frankebel/RAS_DMFT.jl/commit/da721d2607ea8749199a910a9e6d14f8768cc557))

### Fixed

- lognormal broadening on a grid of type `AbstractVector` ([#163](https://github.com/frankebel/RAS_DMFT.jl/pull/163)) ([2cb92dd](https://github.com/frankebel/RAS_DMFT.jl/commit/2cb92dd3b16ae4f3fbd7e28cdb0f578d422ca229))

## [0.10.0] - 2025-12-11

### Changed

- rename package to `RAS_DMFT` ([#161](https://github.com/frankebel/RAS_DMFT.jl/pull/161)) ([fb3e137](https://github.com/frankebel/RAS_DMFT.jl/commit/fb3e137dabb8152441b1fa6cd8466b582e0643dc))

## [0.9.0] - 2025-12-04

### Added

- calculate self-energy using Schur complement ([#131](https://github.com/frankebel/RAS_DMFT.jl/pull/131)) ([caeee45](https://github.com/frankebel/RAS_DMFT.jl/commit/caeee45b179c39f7e412348dd9a15471951b7a8a))
- calculate quasiparticle weight using `quasiparticle_weight` ([#132](https://github.com/frankebel/RAS_DMFT.jl/pull/132)) ([80f7666](https://github.com/frankebel/RAS_DMFT.jl/commit/80f766698407dee54318f4b8a4799cfceae4c035))
- add pole at location 0 with weight 0 ([#137](https://github.com/frankebel/RAS_DMFT.jl/pull/137)) ([c7b6cd8](https://github.com/frankebel/RAS_DMFT.jl/commit/c7b6cd8be714c0d2c43419333fba41fcb4e6e570))
- calculate quasiparticle weight with Gaussian broadening using `quasiparticle_weight_gaussian` ([#145](https://github.com/frankebel/RAS_DMFT.jl/pull/145)) ([ef631cd](https://github.com/frankebel/RAS_DMFT.jl/commit/ef631cde7b0e4174c8220360ba252eed6401a283))
- calculate quasiparticle weight on a grid using `quasiparticle_weight` ([#149](https://github.com/frankebel/RAS_DMFT.jl/pull/149)) ([ba2a2eb](https://github.com/frankebel/RAS_DMFT.jl/commit/ba2a2eb8e074c6723d1b4254c3636edf3ef1d11e))
- merge negative locations for block pole sum ([#150](https://github.com/frankebel/RAS_DMFT.jl/pull/150)) ([133e0ae](https://github.com/frankebel/RAS_DMFT.jl/commit/133e0ae349a3817954a156526f6e0b3a5116964c))
- put block sum of poles on given grid ([#151](https://github.com/frankebel/RAS_DMFT.jl/pull/151)) ([49db03f](https://github.com/frankebel/RAS_DMFT.jl/commit/49db03faeec706372cf94661f574ab2798f23a00))
- flip spectrum for block sum of poles ([#152](https://github.com/frankebel/RAS_DMFT.jl/pull/152)) ([13be8e3](https://github.com/frankebel/RAS_DMFT.jl/commit/13be8e3b5185c90c0f3d9525beba04f76de78519))

### Changed

- `merge_small_poles!` → `merge_small_weight!` ([#130](https://github.com/frankebel/RAS_DMFT.jl/pull/130)) ([35dd38d](https://github.com/frankebel/RAS_DMFT.jl/commit/35dd38d65efaf928a1c0c29cfa4030383c0d8027))
- rename starting wave function ([#133](https://github.com/frankebel/RAS_DMFT.jl/pull/133)) ([709e917](https://github.com/frankebel/RAS_DMFT.jl/commit/709e917b5ca729a6e1a6c71e788f7f7c94486ac6))
  - `starting_Wavefunction` → `Wavefunction_singlet`
  - `starting_CIWavefunction` → `CIWavefunction_singlet`

- `ground_state` → `ground_state!` ([#134](https://github.com/frankebel/RAS_DMFT.jl/pull/134)) ([a34264e](https://github.com/frankebel/RAS_DMFT.jl/commit/a34264eed1918df37d01329d9fb7c3eec335182d))
- require Julia `>=1.12` ([#153](https://github.com/frankebel/RAS_DMFT.jl/pull/153)) ([494ffac](https://github.com/frankebel/RAS_DMFT.jl/commit/494ffac738d7848ec907465dba37a4ad67781bbb))
- logarithmic broadening keeps pole at zero frequency unchanged ([#156](https://github.com/frankebel/RAS_DMFT.jl/pull/156)) ([1b72fba](https://github.com/frankebel/RAS_DMFT.jl/commit/1b72fba8dee742916b357d16f36ef03f90654659))

### Removed

- `self_energy_FG_lorentzian` ([#128](https://github.com/frankebel/RAS_DMFT.jl/pull/128)) ([bc6f5de](https://github.com/frankebel/RAS_DMFT.jl/commit/bc6f5de69c08882285765dcd3d0308f53f855a31))
- `remove_small_poles!` ([#130](https://github.com/frankebel/RAS_DMFT.jl/pull/130)) ([35dd38d](https://github.com/frankebel/RAS_DMFT.jl/commit/35dd38d65efaf928a1c0c29cfa4030383c0d8027))

### Fixed

- addition of `PolesSumBlock` did not create true copies ([#129](https://github.com/frankebel/RAS_DMFT.jl/pull/129)) ([85827bc](https://github.com/frankebel/RAS_DMFT.jl/commit/85827bcbd02fa45563d45fc4c5adbbf591186f62))
- `discretize_similar_weight` did not return sorted poles ([#148](https://github.com/frankebel/RAS_DMFT.jl/pull/148)) ([f1ba260](https://github.com/frankebel/RAS_DMFT.jl/commit/f1ba260d58ea57184e14f8ac4d04a8317397342a))

## [0.8.0] - 2025-07-24

### Added

- `weight` for `Poles` ([#108](https://github.com/frankebel/RAS_DMFT.jl/pull/108)) ([dd5a5da](https://github.com/frankebel/RAS_DMFT.jl/commit/dd5a5da8d0a8704eaa0e5849854a7e0a73713254))
- DMFT with poles on a fixed grid ([#112](https://github.com/frankebel/RAS_DMFT.jl/pull/112)) ([ecf0680](https://github.com/frankebel/RAS_DMFT.jl/commit/ecf06808eb7ec206da841559b7c75ecd2f69c1e2))
- `Core.Array` for `Poles` with complex amplitudes ([#113](https://github.com/frankebel/RAS_DMFT.jl/pull/113)) ([1ba6b19](https://github.com/frankebel/RAS_DMFT.jl/commit/1ba6b19f994693be1dfc6375ca3fdb8cd6ad1985))
- moment of function defined on a grid ([#114](https://github.com/frankebel/RAS_DMFT.jl/pull/114)) ([44c6495](https://github.com/frankebel/RAS_DMFT.jl/commit/44c6495ddb31a1ffe0a7bcf533c66f5dfdd5d862))
- `PolesSum`, `PolesSumBlock` which will replace `Poles` ([#116](https://github.com/frankebel/RAS_DMFT.jl/pull/116)) ([0570d85](https://github.com/frankebel/RAS_DMFT.jl/commit/0570d8567a17704dbdf282ca670b90dbb6b5670f)) ([#117](https://github.com/frankebel/RAS_DMFT.jl/pull/117)) ([4445404](https://github.com/frankebel/RAS_DMFT.jl/commit/44454047e0b661ad23d010ab0e3b2ed225acfcf6))
- `PolesContinuedFraction` ([#118](https://github.com/frankebel/RAS_DMFT.jl/pull/118)) ([7dea995](https://github.com/frankebel/RAS_DMFT.jl/commit/7dea9950644707fbc0211d61075c06197e19a7c4))
- `PolesContinuedFractionBlock` ([#119](https://github.com/frankebel/RAS_DMFT.jl/pull/119)) ([edd8060](https://github.com/frankebel/RAS_DMFT.jl/commit/edd80600f59d27978d055162c85f71dd217fe094))
- conversion between different pole representations ([#123](https://github.com/frankebel/RAS_DMFT.jl/pull/123)) ([f15b907](https://github.com/frankebel/RAS_DMFT.jl/commit/f15b90712b982744dddba506eb68d739389d73bd))

### Changed

- rename `η_gaussian` to `δ_gaussian` ([#110](https://github.com/frankebel/RAS_DMFT.jl/pull/110)) ([c8911ae](https://github.com/frankebel/RAS_DMFT.jl/commit/c8911aefa8a51f1cef2b9f2a143571b3e68c78c6))
- `grid_log` needs positive initial value ([#111](https://github.com/frankebel/RAS_DMFT.jl/pull/111)) ([ebef3fc](https://github.com/frankebel/RAS_DMFT.jl/commit/ebef3fc36df73b4dcd0107c412e724a4966ba8ce))
- `δ_gaussian` works on scalar ([#115](https://github.com/frankebel/RAS_DMFT.jl/pull/115)) ([eed69b5](https://github.com/frankebel/RAS_DMFT.jl/commit/eed69b59c1eae9bebf7d28ad22b123c98c677e48))
- old `Poles` is split into `PolesSum`, `PolesSumBlock` ([#126](https://github.com/frankebel/RAS_DMFT.jl/pull/126)) ([87b2fb5](https://github.com/frankebel/RAS_DMFT.jl/commit/87b2fb5305e0b8f5d396589863889eaf6925b9ca))

### Fixed

- type-stability when taking square root of matrix ([#120](https://github.com/frankebel/RAS_DMFT.jl/pull/120)) ([969d322](https://github.com/frankebel/RAS_DMFT.jl/commit/969d32244b4af02cca57d51b578f958cbd7489a9))
- type-stability for reading in Number ([#121](https://github.com/frankebel/RAS_DMFT.jl/pull/121)) ([1831991](https://github.com/frankebel/RAS_DMFT.jl/commit/18319913a9f8ef5ae01f1c0d782dc4f1dab99f23))

### Removed

- struct `Poles` is split into `PolesSum`, `PolesSumBlock` ([#127](https://github.com/frankebel/RAS_DMFT.jl/pull/127)) ([04189d7](https://github.com/frankebel/RAS_DMFT.jl/commit/04189d786cec43d948c63dfdb5fac8f90e4967a2))

## [0.7.0] - 2025-06-04

### Added

- merge equal poles ([#54](https://github.com/frankebel/RAS_DMFT.jl/pull/54)) ([66b6dc2](https://github.com/frankebel/RAS_DMFT.jl/commit/66b6dc218e681535645e7434da2ae204ffc4bfd3))
- logarithmic grid and log-gaussian broadening ([#55](https://github.com/frankebel/RAS_DMFT.jl/pull/55)) ([f3a09b0](https://github.com/frankebel/RAS_DMFT.jl/commit/f3a09b01e668d3ee188530d1c88ddac1a884b2a6))
- Hubbard-III approximation ([#64](https://github.com/frankebel/RAS_DMFT.jl/pull/64)) ([e5dbed7](https://github.com/frankebel/RAS_DMFT.jl/commit/e5dbed71f02c0949e646a6ba255a7249c2bf8de5))
- warn if wrong spectral weight exists ([#66](https://github.com/frankebel/RAS_DMFT.jl/pull/66)) ([25da3b0](https://github.com/frankebel/RAS_DMFT.jl/commit/25da3b0bf829ba4261094e2cbc1d5bad48bb2924))
- merge poles with small weight to neighbors ([#72](https://github.com/frankebel/RAS_DMFT.jl/pull/72)) ([7ab44a7](https://github.com/frankebel/RAS_DMFT.jl/commit/7ab44a7076250e47c623d3530dc583485d6843c9))
- getters for `Poles` ([#81](https://github.com/frankebel/RAS_DMFT.jl/pull/81)) ([add4ec6](https://github.com/frankebel/RAS_DMFT.jl/commit/add4ec6d421042265f42b8d9745b059a47fa5254))
- weight(s) of each pole ([#82](https://github.com/frankebel/RAS_DMFT.jl/pull/82)) ([b1f1780](https://github.com/frankebel/RAS_DMFT.jl/commit/b1f1780830c4a11b3f340cf3a5a64fc6f40e037b))
- moments of `Poles` ([#83](https://github.com/frankebel/RAS_DMFT.jl/pull/83)) ([1917ce2](https://github.com/frankebel/RAS_DMFT.jl/commit/1917ce246108a72562683d37832b525b1c8e5413))
- similar weight discretization of `Poles` ([#86](https://github.com/frankebel/RAS_DMFT.jl/pull/86)) ([0987324](https://github.com/frankebel/RAS_DMFT.jl/commit/0987324ceed55c9208cc18c6d10633dc12b129e0))
- `Base.issorted` for `Poles` ([#89](https://github.com/frankebel/RAS_DMFT.jl/pull/89)) ([58cd09a](https://github.com/frankebel/RAS_DMFT.jl/commit/58cd09a6f8f8fa19578ca9be45a2c1175d90c42a))
- `Base.allunique` for `Poles` ([#90](https://github.com/frankebel/RAS_DMFT.jl/pull/90)) ([b3d8e52](https://github.com/frankebel/RAS_DMFT.jl/commit/b3d8e52bee8141bc66725e5d5160af5beac5d92d))
- moving poles with negative location to zero ([#91](https://github.com/frankebel/RAS_DMFT.jl/pull/91)) ([5dd9f13](https://github.com/frankebel/RAS_DMFT.jl/commit/5dd9f13db14f88e19f1aa0bae8aadb9e1208666f))
- flip locations of `Poles` ([#95](https://github.com/frankebel/RAS_DMFT.jl/pull/95)) ([5717654](https://github.com/frankebel/RAS_DMFT.jl/commit/5717654fc8d59250b719cd43d4ac9ce99b3d8795)) ([#97](https://github.com/frankebel/RAS_DMFT.jl/pull/97)) ([b560d33](https://github.com/frankebel/RAS_DMFT.jl/commit/b560d330969e40db73357e9a9a727c857c9a5ff1))
- shift locations of `Poles` ([#96](https://github.com/frankebel/RAS_DMFT.jl/pull/96)) ([5b85142](https://github.com/frankebel/RAS_DMFT.jl/commit/5b85142bb4002ed807d1ad730becf295013fe8b3)) ([#97](https://github.com/frankebel/RAS_DMFT.jl/pull/97)) ([0b0ec8f](https://github.com/frankebel/RAS_DMFT.jl/commit/0b0ec8f5f7f075e4ccbd566248887f49ab33c220))
- `correlator_plus()`, `correlator_minus()` ([#98](https://github.com/frankebel/RAS_DMFT.jl/pull/98)) ([3a78427](https://github.com/frankebel/RAS_DMFT.jl/commit/3a78427753859d2a1fc34a9b4af1124be29c1260))
- general block correlator ([#99](https://github.com/frankebel/RAS_DMFT.jl/pull/99)) ([4f20faa](https://github.com/frankebel/RAS_DMFT.jl/commit/4f20faa45b68c969d47491cb859d25ed7c6856a1))

### Changed

- rename `Pole` → `Poles`, `self_energy_pole` → `self_energy_poles` ([#76](https://github.com/frankebel/RAS_DMFT.jl/pull/76)) ([ee98b96](https://github.com/frankebel/RAS_DMFT.jl/commit/ee98b96a051d91be21990a2d2f59300735a798b4))
- rename `move_negative_weight_to_neighbors!` → `merge_negative_weight!` ([#77](https://github.com/frankebel/RAS_DMFT.jl/pull/77)) ([5d88c89](https://github.com/frankebel/RAS_DMFT.jl/commit/5d88c898b2b0507009a78dcacb2f8dac1a36645d))
- rename `merge_equal_poles!` → `merge_degenerate_poles!` ([#79](https://github.com/frankebel/RAS_DMFT.jl/pull/79)) ([c5d2b27](https://github.com/frankebel/RAS_DMFT.jl/commit/c5d2b27ef19d635e0e93e34912c3b0a04668b2ab))
- remove `to_grid_sqr!` ([#80](https://github.com/frankebel/RAS_DMFT.jl/pull/80)) ([cc71c12](https://github.com/frankebel/RAS_DMFT.jl/commit/cc71c12923534a4fe5140930780e737b4a790308))
- calculation of self-energy using `Poles` does not enforce original grid ([#84](https://github.com/frankebel/RAS_DMFT.jl/pull/84)) ([b976fcd](https://github.com/frankebel/RAS_DMFT.jl/commit/b976fcd2b6d1a07b58fe1fb28dc0641d26929e9c))
- calculation of new hybridization function does not enforce original grid ([#85](https://github.com/frankebel/RAS_DMFT.jl/pull/85)) ([bba8809](https://github.com/frankebel/RAS_DMFT.jl/commit/bba8809f196627f18136f2a40582b9383fe15031))
- merge poles that are `tol` apart ([#88](https://github.com/frankebel/RAS_DMFT.jl/pull/88)) ([ad832db](https://github.com/frankebel/RAS_DMFT.jl/commit/ad832dbedb84980f8ea0352af543880ee676b7a3))
- merging negative weight function private ([#93](https://github.com/frankebel/RAS_DMFT.jl/pull/93)) ([5491f81](https://github.com/frankebel/RAS_DMFT.jl/commit/5491f81924d578487557f295cc023f1814e4bcbd))
- remove `update_hybridization_function` for grid ([#102](https://github.com/frankebel/RAS_DMFT.jl/pull/102)) ([dc23810](https://github.com/frankebel/RAS_DMFT.jl/commit/dc23810bb587498d9dff91bfedbec3a459eda58a))
- stabilize Lanczos algorithm ([#103](https://github.com/frankebel/RAS_DMFT.jl/pull/103)) ([607426d](https://github.com/frankebel/RAS_DMFT.jl/commit/607426da8ec4329276bda328bcd33f1c5a171e2d))

### Removed

- DMFT step with Lanczos ([#69](https://github.com/frankebel/RAS_DMFT.jl/pull/69)) ([34633dd](https://github.com/frankebel/RAS_DMFT.jl/commit/34633dd7e88e8e72aad84dbbd496677fd478c434))
- DMFT step with block Lanczos ([#105](https://github.com/frankebel/RAS_DMFT.jl/pull/105)) ([2b7710d](https://github.com/frankebel/RAS_DMFT.jl/commit/2b7710d7c95aaa54ea09e21412f44c20156fa849))
- `solve_impurity` ([#106](http://github.com/frankebel/RAS_DMFT.jl/pull/106)) ([c6c210f](https://github.com/frankebel/RAS_DMFT.jl/commit/c6c210f25e817ddba4c28d9f2b6035e615f86083))

### Fixed

- errors in DMFT self-consistency loop ([#57](https://github.com/frankebel/RAS_DMFT.jl/pull/57)) ([b10573a](https://github.com/frankebel/RAS_DMFT.jl/commit/b10573a5ff1b487878ead85c9fc63c53bc0ed731))
- loss of particle-hole symmetry in natural orbitals ([#58](https://github.com/frankebel/RAS_DMFT.jl/issues/58)) ([#59](https://github.com/frankebel/RAS_DMFT.jl/pull/59)) ([37ad203](https://github.com/frankebel/RAS_DMFT.jl/pull/59/commits/37ad2032a98c06f015ea29152481e9f52333b44c))
- `length` on Pole containing vector and matrix ([#67](https://github.com/frankebel/RAS_DMFT.jl/issues/67)) ([9a6dc41](https://github.com/frankebel/RAS_DMFT.jl/commit/9a6dc418cbeb3d84a074976c3ad15a0fb997513d))

## [0.6.0] - 2025-04-16

### Added

- Green's function and hybridization function on a given grid ([#46](https://github.com/frankebel/RAS_DMFT.jl/pull/46)) ([c30556d](https://github.com/frankebel/RAS_DMFT.jl/commit/c30556ddc816a1b9cf4aa1436d4a3fa88ce6b3fe))
- subtraction of `Pole` ([#47](https://github.com/frankebel/RAS_DMFT.jl/pull/47)) ([8ba021b](https://github.com/frankebel/RAS_DMFT.jl/commit/8ba021bee989cf8ba536fd38a51d43711df62775))
- inversion of `Pole` ([#48](https://github.com/frankebel/RAS_DMFT.jl/pull/48)) ([c34ff69](https://github.com/frankebel/RAS_DMFT.jl/commit/c34ff69349800085824d1430dbd975e83de9e8b3))
- update hybridization function in pole representation ([#50](https://github.com/frankebel/RAS_DMFT.jl/pull/50)) ([c0725f5](https://github.com/frankebel/RAS_DMFT.jl/commit/c0725f5572155657110980cf8caf28ed130a73cb))
- calculate DMFT self-consistency in pure `Pole` representation ([#52](https://github.com/frankebel/RAS_DMFT.jl/pull/52)) ([18c783b](https://github.com/frankebel/RAS_DMFT.jl/commit/18c783bfcada4463ebeb77d85fc8bef7fac357c7))

### Changed

- unify function names of Green's function and hybridization function ([#41](https://github.com/frankebel/RAS_DMFT.jl/pull/41)) ([d095219](https://github.com/frankebel/RAS_DMFT.jl/commit/d095219ade1ae73349ff79e8ea903f69f73159a7)) ([#45](https://github.com/frankebel/RAS_DMFT.jl/pull/45)) ([d4a4bf6](https://github.com/frankebel/RAS_DMFT.jl/commit/d4a4bf6dda05e8a97f9c749bfef08638d1985f89))
- rename functions ([#49](https://github.com/frankebel/RAS_DMFT.jl/pull/49)) ([6fcefff ](https://github.com/frankebel/RAS_DMFT.jl/commit/6fcefffa2f80c817b1dfa95a4001cec880ec6b66))
  - `self_energy` → `self_energy_IFG`
  - `self_energy_gauss` → `self_energy_IFG_gauss`
  - `update_weiss_field` → `update_hybridization_function`

### Fixed

- block Lanczos must return hermitian matrices ([#42](https://github.com/frankebel/RAS_DMFT.jl/pull/42)) ([d0b066a](https://github.com/frankebel/RAS_DMFT.jl/commit/d0b066aba90a8308ea0f9adeece25165e52acaba))
- correlator FR in improved self-energy ([#51](https://github.com/frankebel/RAS_DMFT.jl/pull/51)) ([f374b44](https://github.com/frankebel/RAS_DMFT.jl/commit/f374b444381c9c2a612561cf0d95c32a1733dff3))

## [0.5.1] - 2025-04-02

### Fixed

- documentation ([#40](https://github.com/frankebel/RAS_DMFT.jl/pull/40)) ([032ed29](https://github.com/frankebel/RAS_DMFT.jl/commit/032ed2981c1af41a57eb60616dcbab8f40fc8017))

## [0.5.0] - 2025-04-01

### Added

- calculate Kondo temperature ([#35](https://github.com/frankebel/RAS_DMFT.jl/pull/35)) ([b1fe8ff](https://github.com/frankebel/RAS_DMFT.jl/commit/b1fe8ff94cd895870281b48e2f6a73e0e1c41f7f))
- local Green's function from dispersion relation and optional self-energy ([#36](https://github.com/frankebel/RAS_DMFT.jl/pull/36)) ([925c512](https://github.com/frankebel/RAS_DMFT.jl/commit/925c51201131ee3fc282848aec5be907628ba789))
- partial Green's function ([#37](https://github.com/frankebel/RAS_DMFT.jl/pull/37)) ([1d05c37](https://github.com/frankebel/RAS_DMFT.jl/commit/1d05c37dc7556905ef139266a46ecef003360ecb))
- find chemical potential for desired filling ([#38](https://github.com/frankebel/RAS_DMFT.jl/pull/38)) ([4547335](https://github.com/frankebel/RAS_DMFT.jl/commit/45473351ee71a3a9e736ed836024c35dc97f47ae))
- non-interacting spectral function with Gaussian broadening ([#39](https://github.com/frankebel/RAS_DMFT.jl/pull/39)) ([527b4ab](https://github.com/frankebel/RAS_DMFT.jl/commit/527b4abfcb29952b322d2b1cb569887cf90fa23d))

### Changed

- rename `Greensfunction` to `Pole` to reflect that it is more generic ([#27](https://github.com/frankebel/RAS_DMFT.jl/issues/27)) ([5631a33](https://github.com/frankebel/RAS_DMFT.jl/commit/5631a33405a13b292b0b988edf7b14931b59344a))
- unify IO under `read_hdf5`, `write_hdf5` ([#29](https://github.com/frankebel/RAS_DMFT.jl/pull/29))
- update `Fermions` dependency to `v0.12.0` ([#30](https://github.com/frankebel/RAS_DMFT.jl/pull/30))
- enforce `Pole` pole locations to be on the real axis ([#34](https://github.com/frankebel/RAS_DMFT.jl/pull/34)) ([d094356](https://github.com/frankebel/RAS_DMFT.jl/commit/d094356cf8502aa6a25cdb0049918715182d459f))
- create module `Debug` and clean up namespace ([#34](https://github.com/frankebel/RAS_DMFT.jl/pull/34)) ([c92e6a8](https://github.com/frankebel/RAS_DMFT.jl/commit/c92e6a8a33b5bc787028015a5a40012f7e334985))
