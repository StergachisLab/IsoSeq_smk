BEGIN {
    FS = "\t"
    OFS = "\t"
}

function get_attr(attr, key,    x) {
    if (match(attr, key " \"[^\"]+\"")) {
        x = substr(attr, RSTART, RLENGTH)
        sub(key " \"", "", x)
        sub("\"$", "", x)
        return x
    }
    return ""
}

function flush_block(    i) {
    if (keep_block && n_tx > 0 && n_exon > 0) {
        for (i = 1; i <= n; i++) print block[i]
    }
    n = 0
    n_tx = 0
    n_exon = 0
    keep_block = 0
    current_gid = ""
    delete seen_tx
}

/^#/ || /^[[:space:]]*$/ { next }

NF != 9 { next }

{
    feature = $3
    gid = get_attr($9, "gene_id")
    tid = get_attr($9, "transcript_id")
    gname = get_attr($9, "gene_name")

    if (feature != "gene" && feature != "transcript" && feature != "exon") next
    if ($4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/) next
    if ($7 !~ /^[+-]$/) next

    if (feature == "gene") {
        flush_block()
        if (gid != "" && gname != "") {
            keep_block = 1
            current_gid = gid
            block[++n] = $0
        }
        next
    }

    if (!keep_block || gid != current_gid) next

    if (feature == "transcript") {
        if (tid != "" && gname != "") {
            seen_tx[tid] = 1
            n_tx++
            block[++n] = $0
        }
        next
    }

    if (feature == "exon") {
        if (tid != "" && gname != "" && (tid in seen_tx)) {
            n_exon++
            block[++n] = $0
        }
        next
    }
}

END {
    flush_block()
}
