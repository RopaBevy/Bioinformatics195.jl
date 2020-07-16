module Assignment07

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta

# # uncomment the following line if you intend to use BioSequences types
#using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

function composition(seq)
    seq = normalizeDNA(seq)

    bases = Dict()

    for base in seq
        if haskey(bases, base)
            bases[base] += 1
        else
            bases[base] = 1
        end
    end
    return bases
end

function gc_content(seq)
    bases = composition(seq)

    return (get(bases,'G', 0) + get(bases, 'C', 0)) / length(seq)
end

function complementbase(base)
    base = uppercase(string(base))
    comp = Dict("A"=>'T',
                "T"=>'A',
                "G"=>'C',
                "C"=>'G',
                "N"=>'N')
    return comp[uppercase(base)]
end

function complement(seq)
    seq = normalizeDNA(seq)
    result = ""

    for base in seq
        result = result * complementbase(base)
    end
    return result
end


function reverse_complement(seq)
    seq = normalizeDNA(seq)
    result = ""
    index = length(seq)

    while( index > 0)
        result = result * complementbase(seq[index])
        index = index - 1
    end
    return result
end

function parse_fasta(path)
    heading = Vector()
    seqs = Vector()
    sequence = ""

    for line in eachline(path)
        if startswith(line, '>') 
            if length(sequence) != 0
                push!(seqs, sequence)
            end
            sequence = "" 
            push!(heading, line[2:end])
        else 
            line = normalizeDNA(line)
            sequence = sequence * line
        end
    end
    if length(sequence) != 0
        push!(seqs, sequence)
    end

    return heading, seqs
end

end # module Assignment07