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
    seq = normalizeDNA(seq)
    count = 0

    for base in seq
        if base == DNA_G || base == DNA_C
            count = count + 1
        end
    end
    return count / length(seq)
end

function complement(seq)
    seq = normalizeDNA(seq)
    result = ""

    for base in seq
        if base == "A"
            result = result * "T"
        elseif base == "T"
            result = result * "A"
        elseif base == "G"
            result = result * "C"
        elseif base == "C"
            result = result * "G"
        else result = result * "N"
        end 
    end
    return result
end

function reverse_complement(seq)
    seq = normalizeDNA(seq)
    result = ""
    index = length(seq)

    while( index > 0)
        if seq[index] == "A"
            result = result * "T"
        elseif seq[index] == "T"
            result = result * "A"
        elseif seq[index] == "G"
            result = result * "C"
        elseif seq[index] == "C"
            result = result * "G"
        else result = result * "N"
        end 
        index = index - 1
    end
    return result
end

function parse_fasta(path)
    vector1 = Vector()
    vector2 = Vector()
    sequence = ""

    for line in eachline(path)
        if startswith(line, '>') 
            if length(sequence) != 0
                push!(vector2, sequence)
            end
            sequence = "" 
            push!(vector1, line[2:end])
        else sequence = sequence * line
        end
    end
    if length(sequence) != 0
        push!(vector2, sequence)
    end

    return vector1, vector2
end

end # module Assignment07