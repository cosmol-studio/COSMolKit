# Third-party notices for `cosmolkit-fingerprints`

This file records upstream source reproduced or quoted by this crate. It
supplements the COSMolKit project license; it does not replace or change that
license, relicense the project, or change package metadata.

## RDKit

RDKit source pin: `351f8f378f8ad6bbd517980c38896e66bf907af8c`.

The following local modules contain source-derived code or source anchors:

| Local module | Pinned RDKit source/header pointers |
|---|---|
| `src/hash.rs` | `Code/RDGeneral/hash/hash_fwd.hpp`, `hash.hpp`, `extensions.hpp` |
| `src/morgan.rs` | `Code/GraphMol/Fingerprints/MorganFingerprints.cpp`, `MorganGenerator.cpp`, `MorganGenerator.h`, `FingerprintGenerator.cpp`, `FingerprintGenerator.h`; `Code/RDGeneral/RDProps.h`, `Dict.h`, `RDValue.h`, and `hash/{hash.hpp,hash_fwd.hpp,extensions.hpp}` |
| `src/generator.rs` | `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp`, `FingerprintGenerator.h`, `MorganGenerator.cpp`, `MorganGenerator.h`, and `MorganFingerprints.cpp` |
| `src/additional_output.rs` | `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp`, `FingerprintGenerator.h`, and `Wrap/FingerprintGeneratorWrapper.cpp` |
| `src/prepared.rs` | `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp`; `Code/GraphMol/ROMol.cpp`, `ROMol.h`, `MolOps.cpp`, `MolOps.h`, `Chirality.cpp`, and `Chirality.h` |
| `src/invariants.rs` | `Code/GraphMol/Fingerprints/FingerprintUtil.cpp`, `FingerprintUtil.h`, `MorganGenerator.cpp`, and `MorganGenerator.h`; `Code/GraphMol/Atom.cpp`, `Atom.h`, `ROMol.cpp`, `ROMol.h`, `PeriodicTable.cpp`, `PeriodicTable.h`, `RingInfo.cpp`, and `RingInfo.h`; `Code/RDGeneral/RDProps.h`, `Dict.h`, `RDValue.h`, and `hash/{hash.hpp,hash_fwd.hpp,extensions.hpp}` |

Copyright notices in the copied RDKit source closure include:

- `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp` and `.h`:
  Copyright (C) 2018-2025 Boran Adas and other RDKit contributors.
- `Code/GraphMol/Fingerprints/MorganGenerator.cpp`:
  Copyright (C) 2018-2025 Boran Adas and other RDKit contributors.
- `Code/GraphMol/Fingerprints/MorganGenerator.h`:
  Copyright (C) 2018-2022 Boran Adas and other RDKit contributors.
- `Code/GraphMol/Fingerprints/FingerprintUtil.cpp`:
  Copyright (C) 2018-2025 Boran Adas and other RDKit contributors.
- `Code/GraphMol/Fingerprints/FingerprintUtil.h`:
  Copyright (C) 2018 Boran Adas, Google Summer of Code.
- `Code/GraphMol/Atom.cpp` and `.h`: Copyright (C) 2001-2024 Greg
  Landrum and other RDKit contributors.
- `Code/GraphMol/ROMol.cpp`: Copyright (C) 2003-2024 Greg Landrum and other
  RDKit contributors; `ROMol.h`: Copyright (C) 2003-2022 Greg Landrum and
  other RDKit contributors.
- `Code/GraphMol/MolOps.cpp`: Copyright (C) 2001-2023 Greg Landrum and other
  RDKit contributors; `MolOps.h`: Copyright (C) 2001-2024 Greg Landrum and
  other RDKit contributors.
- `Code/GraphMol/Chirality.cpp`: Copyright (C) 2004-2024 Greg Landrum and
  other RDKit contributors; `Chirality.h`: Copyright (C) 2008-2022 Greg
  Landrum and other RDKit contributors.
- `Code/GraphMol/PeriodicTable.cpp`: Copyright (C) 2001-2006 Rational
  Discovery LLC; `PeriodicTable.h`: Copyright (C) 2001-2011 Rational
  Discovery LLC.
- `Code/GraphMol/RingInfo.cpp`: Copyright (C) 2004-2019 Greg Landrum and
  Rational Discovery LLC; `RingInfo.h`: Copyright (C) 2004-2022 Greg Landrum
  and other RDKit contributors.
- `Code/RDGeneral/RDProps.h`: Copyright (C) 2016-2026 Brian Kelley and other
  RDKit contributors.
- `Code/RDGeneral/Dict.h`: Copyright (C) 2003-2026 Greg Landrum and other
  RDKit contributors.
- `Code/RDGeneral/hash/hash.hpp`, `hash_fwd.hpp`, and `extensions.hpp`:
  Copyright 2005-2009 Daniel James; `hash.hpp` also records modifications by
  Greg Landrum to obtain portable hashes across machines.

### RDKit BSD 3-Clause License

The pinned source's root `license.txt` contains this complete notice:

```text
BSD 3-Clause License

Copyright (c) 2006-2015, Rational Discovery LLC, Greg Landrum, and Julie Penzotti and others
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### `MorganFingerprints.cpp` notice

The pinned file carries this separate complete notice:

```text
Copyright (c) 2009-2022, Novartis Institutes for BioMedical Research Inc.
and other RDKit contributors
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the name of Novartis Institutes for BioMedical Research Inc.
      nor the names of its contributors may be used to endorse or promote
      products derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### `RDValue.h` notice

The pinned header carries this additional copyright and complete BSD notice:

```text
Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the name of Novartis Institutes for BioMedical Research Inc.
      nor the names of its contributors may be used to endorse or promote
      products derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Boost.Random

Boost.Random source pin: `09594050246adcf944be2b2c8dd2b9901083de79`.
The exact RNG comments in `src/rng.rs` are derived from these immutable
headers:

| Upstream header | Header copyright notice |
|---|---|
| `include/boost/random/mersenne_twister.hpp` | Jens Maurer 2000-2001; Steven Watanabe 2010 |
| `include/boost/random/uniform_int.hpp` | Jens Maurer 2000-2001 |
| `include/boost/random/uniform_int_distribution.hpp` | Jens Maurer 2000-2001; Steven Watanabe 2011 |
| `include/boost/random/variate_generator.hpp` | Jens Maurer 2002; Steven Watanabe 2011 |
| `include/boost/random/detail/seed.hpp` | Steven Watanabe 2009 |
| `include/boost/random/detail/signed_unsigned_tools.hpp` | Jens Maurer 2006 |
| `include/boost/random/detail/ptr_helper.hpp` | Jens Maurer 2002 |

## Boost.DynamicBitset

Boost.DynamicBitset source pin: `8e20aa1462bf6dcadc338835df529a6d568431b1`.
`src/morgan.rs` reproduces selected dynamic-bitset operations from
`include/boost/dynamic_bitset/dynamic_bitset.hpp`, whose notices identify
Chuck Allison and Jeremy Siek (2001-2002), Gennaro Prota (2003-2006, 2008),
Ahmed Charles (2014), Glen Joseph Fernandes (2014), Riccardo Marcangelo
(2014), and Evgeny Shulgin (2018).

## Boost Software License 1.0

The following is the complete `LICENSE_1_0.txt` text from the Boost 1.85.0
source distribution. It applies to the Boost.Random and Boost.DynamicBitset
source-derived portions above, and to the RDKit modified hash headers where
those headers carry the same license.

```text
Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
```

## Avalon Toolkit

`src/avalon/` reproduces the pinned AvalonToolkit_2.0.5-pre.3 parser and
fingerprint engine (RDKit ava-formake distribution, archive MD5
`7a20c25a7e79f3344e0f9f49afa03351`). Engine source files also carry:
Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
All rights reserved.

The complete distribution license follows:

```text
Copyright 2001-2011 Novartis Pharma AG. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

   3. Neither the name of Novartis Institutes for BioMedical Research Inc.
      nor the names of its contributors may be used to endorse or promote
      products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY NOVARTIS PHARMA AG ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of Novartis Pharma AG.
```
