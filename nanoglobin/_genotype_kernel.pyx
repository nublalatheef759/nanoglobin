# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
"""Required Cython kernel for assay-observable genotype scoring.

Read/product edit distance already runs in RapidFuzz's native implementation. This
kernel handles the remaining class-by-evidence likelihood, product-composition and
required-product terms. Equivalent molecule observations are grouped before entry.
"""

import numpy as np
cimport numpy as cnp

from libc.math cimport lgamma, log

cnp.import_array()


def score_dense_classes(
    cnp.ndarray[cnp.uint8_t, ndim=2] class_exact,
    cnp.ndarray[cnp.uint8_t, ndim=2] class_product,
    cnp.ndarray[cnp.uint8_t, ndim=2] class_required,
    cnp.ndarray[cnp.int64_t, ndim=2] class_multiplicity,
    cnp.ndarray[cnp.uint64_t, ndim=2] class_haplotype_masks,
    cnp.ndarray[cnp.int64_t, ndim=1] group_product,
    cnp.ndarray[cnp.int64_t, ndim=1] group_key,
    cnp.ndarray[cnp.uint64_t, ndim=2] group_candidate_masks,
    cnp.ndarray[cnp.float64_t, ndim=1] group_weight,
    cnp.ndarray[cnp.float64_t, ndim=1] group_edit_penalty,
    cnp.ndarray[cnp.float64_t, ndim=1] observed_counts,
    cnp.ndarray[cnp.float64_t, ndim=1] product_efficiencies,
    cnp.ndarray[cnp.float64_t, ndim=1] log_priors,
    double artifact_probability,
    double product_background_mass,
    double count_concentration,
    double minimum_dirichlet_alpha,
    double dropout_probability,
):
    """Score every observable genotype class using typed dense arrays."""

    cdef Py_ssize_t class_count = class_exact.shape[0]
    cdef Py_ssize_t key_count = class_exact.shape[1]
    cdef Py_ssize_t product_count = class_product.shape[1]
    cdef Py_ssize_t group_count = group_product.shape[0]
    cdef Py_ssize_t mask_word_count = class_haplotype_masks.shape[1]
    cdef Py_ssize_t c, g, p, word
    cdef cnp.int64_t product_index, key_index
    cdef double read_value, count_value, dropout_value, total_value
    cdef double log_probability, edit_penalty, mass, total_mass
    cdef double alpha, alpha0, n_total, count
    cdef double log_exact = log(1.0 - artifact_probability)
    cdef double log_half = log(
        artifact_probability
        if artifact_probability > 0.5 * (1.0 - artifact_probability)
        else 0.5 * (1.0 - artifact_probability)
    )
    cdef double log_quarter = log(
        artifact_probability
        if artifact_probability > 0.25 * (1.0 - artifact_probability)
        else 0.25 * (1.0 - artifact_probability)
    )
    cdef double log_artifact = log(artifact_probability)
    cdef double log_dropout = log(dropout_probability)
    cdef double log_present = log(1.0 - dropout_probability)
    cdef bint product_present, exact_present, candidate_overlap

    if class_product.shape[0] != class_count:
        raise ValueError("class_product row count does not match class_exact")
    if class_required.shape[0] != class_count:
        raise ValueError("class_required row count does not match class_exact")
    if class_multiplicity.shape[0] != class_count:
        raise ValueError("class_multiplicity row count does not match class_exact")
    if class_required.shape[1] != product_count:
        raise ValueError("class_required product count does not match")
    if class_multiplicity.shape[1] != product_count:
        raise ValueError("class_multiplicity product count does not match")
    if class_haplotype_masks.shape[0] != class_count or mask_word_count < 1:
        raise ValueError("class_haplotype_masks shape does not match")
    if observed_counts.shape[0] != product_count:
        raise ValueError("observed_counts length does not match product count")
    if product_efficiencies.shape[0] != product_count:
        raise ValueError("product_efficiencies length does not match product count")
    if log_priors.shape[0] != class_count:
        raise ValueError("log_priors length does not match class count")
    if group_key.shape[0] != group_count:
        raise ValueError("group_key length does not match")
    if group_candidate_masks.shape[0] != group_count:
        raise ValueError("group_candidate_masks row count does not match")
    if group_candidate_masks.shape[1] != mask_word_count:
        raise ValueError("candidate and class haplotype masks do not match")
    if group_weight.shape[0] != group_count:
        raise ValueError("group_weight length does not match")
    if group_edit_penalty.shape[0] != group_count:
        raise ValueError("group_edit_penalty length does not match")

    cdef cnp.ndarray[cnp.float64_t, ndim=1] read_out = np.zeros(
        class_count, dtype=np.float64
    )
    cdef cnp.ndarray[cnp.float64_t, ndim=1] count_out = np.zeros(
        class_count, dtype=np.float64
    )
    cdef cnp.ndarray[cnp.float64_t, ndim=1] dropout_out = np.zeros(
        class_count, dtype=np.float64
    )
    cdef cnp.ndarray[cnp.float64_t, ndim=1] score_out = np.zeros(
        class_count, dtype=np.float64
    )

    with nogil:
        for c in range(class_count):
            read_value = 0.0
            for g in range(group_count):
                product_index = group_product[g]
                key_index = group_key[g]
                product_present = (
                    product_index >= 0
                    and product_index < product_count
                    and class_product[c, product_index] != 0
                )
                exact_present = (
                    key_index >= 0
                    and key_index < key_count
                    and class_exact[c, key_index] != 0
                )
                candidate_overlap = False
                for word in range(mask_word_count):
                    if (
                        class_haplotype_masks[c, word]
                        & group_candidate_masks[g, word]
                    ) != 0:
                        candidate_overlap = True
                        break
                if exact_present:
                    log_probability = log_exact
                elif key_index == -1 and product_present and candidate_overlap:
                    # -1 means no compiled-sequence hash. -2 means a supplied
                    # hash is absent from every candidate class.
                    log_probability = log_half
                elif product_present:
                    log_probability = log_quarter
                else:
                    log_probability = log_artifact
                edit_penalty = group_edit_penalty[g]
                read_value += group_weight[g] * log_probability - edit_penalty

            total_mass = 0.0
            alpha0 = 0.0
            n_total = 0.0
            for p in range(product_count):
                if observed_counts[p] > 0.0 or class_multiplicity[c, p] > 0:
                    mass = (
                        class_multiplicity[c, p] * product_efficiencies[p]
                        + product_background_mass
                    )
                    total_mass += mass
            count_value = 0.0
            if total_mass > 0.0:
                for p in range(product_count):
                    if observed_counts[p] > 0.0 or class_multiplicity[c, p] > 0:
                        mass = (
                            class_multiplicity[c, p] * product_efficiencies[p]
                            + product_background_mass
                        )
                        alpha = count_concentration * mass / total_mass
                        if alpha < minimum_dirichlet_alpha:
                            alpha = minimum_dirichlet_alpha
                        alpha0 += alpha
                        n_total += observed_counts[p]
                count_value = lgamma(alpha0) - lgamma(alpha0 + n_total)
                for p in range(product_count):
                    if observed_counts[p] > 0.0 or class_multiplicity[c, p] > 0:
                        mass = (
                            class_multiplicity[c, p] * product_efficiencies[p]
                            + product_background_mass
                        )
                        alpha = count_concentration * mass / total_mass
                        if alpha < minimum_dirichlet_alpha:
                            alpha = minimum_dirichlet_alpha
                        count = observed_counts[p]
                        count_value += lgamma(alpha + count) - lgamma(alpha)

            dropout_value = 0.0
            for p in range(product_count):
                if class_required[c, p] != 0:
                    if observed_counts[p] > 0.0:
                        dropout_value += log_present
                    else:
                        dropout_value += log_dropout

            total_value = read_value + count_value + dropout_value + log_priors[c]
            read_out[c] = read_value
            count_out[c] = count_value
            dropout_out[c] = dropout_value
            score_out[c] = total_value

    return read_out, count_out, dropout_out, score_out
