"""Decode ConfSeq strings with the reference and fast template backends."""

from __future__ import annotations

import cosmolkit


# Tokenized corpus record: dude_aa2ar__L021087 candidate 0.
# Spaces are ConfSeq token boundaries and must be preserved.
confseq = (
    "N C <112> | ( = O ) <84> C <111> | <-45> n 1 c c ( <21> N <123> | "
    "<173> C <116> | ( = O ) <120> c 2 c c c c c 2 <112> C <112> | <-66> "
    "N 2 <-172> C ( = O ) <-2> C <0> N <3> C <174> 2 = O ) c n 1"
)

reference = cosmolkit.confseq.decode(
    confseq,
    optimize_with_uff=False,
    template_backend="distance_geometry",
)

fast = cosmolkit.confseq.decode(
    confseq,
    optimize_with_uff=False,
    template_backend="fast_geometry",
)

print("Corpus record: dude_aa2ar__L021087 candidate 0")
print("ConfSeq numeric tokens:", confseq.count("<"))
print("DistanceGeometry decoded atoms:", reference.num_atoms(), "conformers:", reference.num_conformers())
print("FastGeometry decoded atoms:", fast.num_atoms(), "conformers:", fast.num_conformers())
