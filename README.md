LineageResolver
Joint probabilistic adjudication of ambiguous cell identity in single-cell RNA-seq


LineageResolver is a tool for resolving ambiguous immune cell identities in single-cell RNA sequencing data.

Many immune cell populations share highly similar transcriptional programs, making them difficult to distinguish using gene expression alone. For example, cytotoxic NK cells, γδ T cells, and cytotoxic CD8 T cells often cluster together in transcriptomic space.

Current analysis workflows treat ambient RNA correction, RNA-based annotation, and orthogonal evidence (such as TCR fragments or protein markers) as separate steps. This can lead to confident but incorrect cell type assignments.

LineageResolver performs joint probabilistic inference over cell identity by modeling RNA expression, ambient contamination, and partial orthogonal evidence simultaneously.
