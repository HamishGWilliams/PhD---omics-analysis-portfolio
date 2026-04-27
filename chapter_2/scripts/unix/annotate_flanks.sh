awk -v upstream=1000 -v downstream=1000 '
BEGIN { FS=OFS="\t" }

function get_attr(attrs, key,   n,i,f,pair) {
    n = split(attrs, f, ";")
    for (i = 1; i <= n; i++) {
        split(f[i], pair, "=")
        if (pair[1] == key) return pair[2]
    }
    return ""
}

{
    if ($0 ~ /^#/) {
        print
        next
    }

    if ($3 == "gene") {
        gene_id = get_attr($9, "ID")
        chr     = $1
        start   = $4
        end     = $5
        strand  = $7

        if (gene_id != "") {

            # promoter region
            if (strand == "+") {
                p_start = start - upstream
                if (p_start < 1) p_start = 1
                p_end = start - 1
            } else if (strand == "-") {
                p_start = end + 1
                p_end   = end + upstream
            }

            if (p_start <= p_end) {
                printf "%s\tcustom\tpromoter\t%d\t%d\t.\t%s\t.\tID=promoter_of_%s;Parent=%s\n", chr, p_start, p_end, strand, gene_id, gene_id
            }

            # downstream region based only on gene boundary
            if (strand == "+") {
                d_start = end + 1
                d_end   = end + downstream
            } else if (strand == "-") {
                d_start = start - downstream
                if (d_start < 1) d_start = 1
                d_end   = start - 1
            }

            if (d_start <= d_end) {
                printf "%s\tcustom\tdownstream_region\t%d\t%d\t.\t%s\t.\tID=downstream_of_%s;Parent=%s\n", chr, d_start, d_end, strand, gene_id, gene_id
            }
        }
    }

    print
}
' A_Equina.gff3 > A_Equina.gff3 