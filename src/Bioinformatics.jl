module Bioinformatics

greet() = print("Hello World!")

complement(nucleotide) = Dict('A'=>'T','T'=>'A','G'=>'C','C'=>'G')[nucleotide]


function reversecomplement(pattern)
    r = Char[]
    for i=1:length(pattern)
        nucleotide = pattern[end-i+1]
        push!(r, complement(nucleotide))
    end
    return string(r...)
end


function patterncount(text, pattern)
    count = 0
    for i = 1:(length(text)-length(pattern))+1
        if text[i:i+length(pattern)-1] == pattern
            count += 1
        end
    end
    return count
end


function frequentwords(text, k)
    count = zeros(length(text))
    frequentpatterns = Set()
    for i = 1:length(text)-k+1
        pattern = String(text[i:i+k-1])
        count[i] = patterncount(text, pattern)
    end
    maxCount = maximum(count)
    for i = 1:length(text)-k+1
        if count[i] == maxCount
            push!(frequentpatterns, text[i:i+k-1])
        end
    end
    return frequentpatterns
end


function patterntonumber(pattern)
    d=Dict('A'=>"00",'C'=>"01",'G'=>"10",'T'=>"11")
    binrep=string([d[bp] for bp in pattern]...)
    parse(Int,binrep,base=2)
end


function patterntonumber2(pattern)
    if length(pattern)==1
        return Dict("A"=>0,"C"=>1,"G"=>2,"T"=>3)[pattern]
    end
    return 4*patterntonumber2(pattern[1:end-1]) + patterntonumber2(pattern[end:end])
end


function patterntonumber(pattern)
    d=Dict('A'=>"00",'C'=>"01",'G'=>"10",'T'=>"11")
    binrep=string([d[bp] for bp in pattern]...)
    parse(Int,binrep,base=2)
end


function numbertopattern(number, k)
    r = String[]
    d = Dict(0=>"A",1=>"C",2=>"G",3=>"T")
    while number >= 1
        _rm = rem(number, 4)
        number = div(number, 4)
        push!(r, d[_rm])
    end
    while length(r) < k
        push!(r, "A")
    end
    return string(reverse(r)...)
end


function computingfrequencies(text, k)
    frequencyarray=zeros(4^k)
    for i = 1:length(text)-k+1
        pattern = text[i:i+k-1]
        j = patterntonumber(pattern) + 1
        frequencyarray[j] = frequencyarray[j] + 1
    end
    return frequencyarray
end


function fasterfrequentwords(text, k)
    frequentpatterns = Set()
    frequencyarray = computingfrequencies(text, k)
    maxCount = maximum(frequencyarray)
    for i=1:4^k
        if frequencyarray[i] == maxCount
            pattern = numbertopattern(i, k)
            push!(frequentpatterns, pattern)
        end
    end
    return frequentpatterns
end


function findingfrequentwordsbysorting(text, k)
    frequentpatterns = Set()
    count = zeros(length(text))
    index = zeros(length(text))
    for i = 1:length(text)-k+1
        pattern = text[i:i+k-1]
        index[i] = patterntonumber(pattern)
        count[i] = 1
    end
    sortedindex = sort(index)
    for i = 2:length(text)-k+1
        if sortedindex[i] == sortedindex[i-1]
            count[i] = count[i-1] + 1
        end
    end
    maxCount = maximum(count)
    for i = 1:length(text)-k+1
        if count[i] == maxCount
            pattern = numbertopattern(sortedindex[i], k)
            push!(frequentpatterns,pattern)
        end
    end
    return frequentpatterns
end


function clumpfinding(genome, k, l, t)
    frequentpatterns = Set()
    clump = zeros(4^k)
    for i = 1:length(genome)-l
        text = genome[i:i+l-1]
        frequencyarray = computingfrequencies(text, k)
        for index = 1:4^k
            if frequencyarray[index] >= t
                clump[index] = 1
            end
        end
    end
    for i = 1:4^k
        if clump[i] == 1
            pattern = numbertopattern(i-1, k)
            push!(frequentpatterns,pattern)
        end
    end
    return frequentpatterns
end


function betterclumpfinding(genome, k, l, t)
    frequentpatterns = Set()
    clump = zeros(4^k)
    text = genome[1:l]
    frequencyarray = computingfrequencies(text, k)
    for i=1:4^k
        if frequencyarray[i] >= t
            clump[i] = 1
        end
    end
    for i=2:length(genome)-l
        firstpattern = genome[i-1:i+k-2]
        index = patterntonumber(firstpattern) + 1
        frequencyarray[index] -= 1
        lastpattern = genome[i+l-k:i+l-1]
        index = patterntonumber(lastpattern) + 1
        frequencyarray[index] += 1
        if frequencyarray[index] >= t
            clump[index] = 1
        end
    end
    for i=1:4^k
        if clump[i] == 1
            pattern = numbertopattern(i-1, k)
            push!(frequentpatterns, pattern)
        end
    end
    return frequentpatterns
end


function hamming(str1,str2)
    d = 0
    for i=1:length(str1)
        if str1[i]!=str2[i]
           d += 1
        end
    end
    return d
end


function skew(genome)
    r = 0
    s = zeros(length(genome)+1)
    for i=1:length(genome)
        if genome[i]=='C'
            r -=1
        elseif genome[i]=='G'
            r +=1
        end
        s[i+1] = r
    end
    return s
end


function minimumskew(genome)
    s = skew(genome)
    minindices = Array{Int,1}()
    for i = 1:length(genome)
        if s[i] == minimum(s)
            push!(minindices, i)
        end
    end
    return minindices
end


function approximatepatternmatch(pattern, text, d)
    matchindices = Array{Int,1}()
    for i = 1:(length(text)-length(pattern))+1
        if hamming(text[i:i+length(pattern)-1],pattern)<=d
            push!(matchindices, i)
        end
    end
    return matchindices
end


function approximatepatterncount(text, pattern, d)
    count = 0
    for i = 1:(length(text)-length(pattern))+1
        if hamming(text[i:i+length(pattern)-1],pattern)<=d
            count += 1
        end
    end
    return count
end


function immediateneighbours(pattern)
    neighbourhood = Set([pattern])
    for i=1:length(pattern)
        sym = pattern[i]
        for nuc in ['A','T','G','C']
            if nuc == sym
                continue
            end
            neighbour = string(pattern[1:i-1],nuc,pattern[i+1:end])
            push!(neighbourhood, neighbour)
        end
    end
    return neighbourhood
end


function neighbours(pattern, d)
    if d == 0
        return [pattern]
    end
    if length(pattern) == 1
        return Set(["A", "C", "G", "T"])
    end
    neighbourhood = Set{String}()
    suffixneighbours = neighbours(pattern[2:end], d)
    for text in suffixneighbours
        if hamming(pattern[2:end],text) < d
            for nuc in ['A', 'C', 'G', 'T']
                push!(neighbourhood, string(nuc,text))
            end
        else
            push!(neighbourhood, string(pattern[1], text))
        end
    end
    return neighbourhood
end


function frequentwordswithmismatches(text, k, d)
    frequentpatterns = Set{String}()
    neighbourhoods = String[]
    for i = 1:length(text)-k+1
        for nb in neighbours(text[i:i+k-1], d)
            push!(neighbourhoods, nb)
        end
    end
    count = ones(length(neighbourhoods))
    index = zeros(length(neighbourhoods))
    for i = 1:length(neighbourhoods)
        pattern = neighbourhoods[i] 
        index[i] = patterntonumber(pattern)
    end
    sortedindex = sort(index)
    for i = 1:length(neighbourhoods) - 1
        if sortedindex[i] == sortedindex[i + 1]
            count[i + 1] = count[i] + 1
        end
    end
    maxcount = maximum(count)
    for i = 1:length(neighbourhoods)
        if count[i] == maxcount
            pattern = numbertopattern(sortedindex[i], k)
            push!(frequentpatterns, pattern)
        end
    end
    return frequentpatterns
end


function motifenumeration(k, d, dna)
    patterns = Set{String}()
    for i=1:length(dna[1])-k+1
        pattern = dna[1][i:i+k-1]
        for pattern2 in neighbours(pattern, d)
            global iscommon = true
            for j=1:length(dna)
                if approximatepatterncount(dna[j], pattern2, d) < 1
                    global iscommon = false
                end
            end
            if iscommon
                push!(patterns, pattern2)
            end
        end
    end
    return patterns
end

"""
STOP and Think: What are the maximum and minimum possible values for the entropy of a probability distribution containing four values?
maximum entropy = sum_4_options - 0.25 log2 0.25 = log2 4 = 2 bits
minimum entropy = - 0.25 * log2 1 = 0.25 * 0 = 0 bits
"""
motifs = 
['T' 'C' 'G' 'G' 'G' 'G' 'G' 'T' 'T' 'T' 'T' 'T'
 'C' 'C' 'G' 'G' 'T' 'G' 'A' 'C' 'T' 'T' 'A' 'C'
 'A' 'C' 'G' 'G' 'G' 'G' 'A' 'T' 'T' 'T' 'T' 'C'
 'T' 'T' 'G' 'G' 'G' 'G' 'A' 'C' 'T' 'T' 'T' 'T'
 'A' 'A' 'G' 'G' 'G' 'G' 'A' 'C' 'T' 'T' 'C' 'C'
 'T' 'T' 'G' 'G' 'G' 'G' 'A' 'C' 'T' 'T' 'C' 'C'
 'T' 'C' 'G' 'G' 'G' 'G' 'A' 'T' 'T' 'C' 'A' 'T'
 'T' 'C' 'G' 'G' 'G' 'G' 'A' 'T' 'T' 'C' 'C' 'T'
 'T' 'A' 'G' 'G' 'G' 'G' 'A' 'A' 'C' 'T' 'A' 'C'
 'T' 'C' 'G' 'G' 'G' 'T' 'A' 'T' 'A' 'A' 'C' 'C']
p = vcat([sum([x==i for x in motifs],dims=1) for i in ['A','C','G','T']]...)./size(motifs,1)
h = -p .* log2.(p)
h[isnan.(h)] .= 0
h = sum(h,dims=1)
H = sum(h)

function distancebetweenpatternsandstrings(pattern, dna)
    k = length(pattern)
    distance = 0
    for string in dna
        hammingdistance = Inf
        for i=1:length(string)-k+1
            pattern2 = string[i:i+k-1]
            d = hamming(pattern, pattern2)
            if hammingdistance > d
                hammingdistance = d
            end
        end
        distance += hammingdistance
    end
    return distance
end

function medianstring(dna, k)
    distance = Inf
    median = 'A'^k
    for i=1:4^k
        pattern = numbertopattern(i, k)
        d = distancebetweenpatternsandstrings(pattern, dna)
        if distance > d
            distance = d
            median = pattern
        end
    end
    return median
end

function profileprobability(pattern, profile)
    p = 1
    row = Dict('A'=>1,'C'=>2,'G'=>3,'T'=>4)
    for (j,i) in enumerate([row[x] for x in pattern])
        p *= profile[i,j]
    end
    return p
end

function profilemostprobable(text, k, profile)
    bestp = 0
    besti = 1
    for i=1:length(text)-k+1
        pattern = text[i:i+k-1]
        p = profileprobability(pattern, profile)
        if p > bestp
            bestp = p
            besti = i
        end
    end
    bestpattern = text[besti:besti+k-1]
end

text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
k = 5
profile = [0.2 0.2 0.3 0.2 0.3
 0.4 0.3 0.1 0.5 0.1
 0.3 0.3 0.5 0.2 0.4
 0.1 0.2 0.1 0.1 0.2]

text = "CGCGATGTGGGGTCTGCCACTCAGCTGGTAAACGCAAACCACGCAGCTGTCCTGCAGGACACGTCGATATGAAGCTATAAGTGACATGGCAGGCACTATACTGTCCCCTTCCTTGTCCAATCGATGACGGTGACGCACAATCTTATGTTCTGCATCGGCTCACTCGGCGTCGGATGGCTCGCCGACAACGAGACTTCCGTCAAGTCGATGAGGAAAAATAGTCTCGGCGTACAGATCTAAACCCACAAGTCCGTTTGGTATACGCAGTCCACCACGTCCATAGACGACCTCGCAATTTAAACCCTTCAGGACAGACAGTAGGGATATTAGAGTTCCGGCGACAAGATCCTCTGAGGTGTCGATTCGCCACCGATTGCAACTGCTACCTACCGACGGCTCTAGGATCGTAGATTTGAGAAACCGATACGCGGGAGCAAGCCTATGGAAAGACTTACTTCATCTTAGAGGGATTATAGAAGCAACTTTATAGTCTCGCCACCGAATGTAATCGTTAGCATGGTGGTACGGCTCGCTGTGTCTGCCAACCCTAAACTGGTTATATGTTAGGCACACGTAACTTGAGGATAATTATACAAATGAATCCATCCGTCTCGCTGTCTCGACGGCTGGAGATGTATCCAAACTTAGAAAACAGAGATTGAGTTTCCAGAACCTCAAGCCCGAGCGCCACATTGATTTCAATTCCGTTTGCAGAGGACATCTGAAAGATGCAATCCGCTCCTGAATGCGCGATAGTCTAAAATCTGGAATCAGAAAGTCCTTCTCCAGAGGATCGGCAAGCTTTTCCTAGACCTCCGGCCTCCGGCGGTCAGAAACAGAAGCCTTGACCGGACGTGAAAATCCCACAAAAAAAGGACCTGATTTACTCTAGATCTTTGCGTCACCTCTCTAGTAGGCCGAACAACCACCTTTAATCGACTCATGGACTATTATCTGCGGAATCTGAACAACACTCCAGCCATTCCGCATCAATAGCTTC"
k = 15
profile = [
0.303 0.273 0.212 0.212 0.212 0.273 0.318 0.242 0.227 0.303 0.242 0.242 0.242 0.273 0.242
0.212 0.197 0.303 0.288 0.227 0.182 0.212 0.106 0.242 0.197 0.288 0.288 0.288 0.152 0.197
0.182 0.288 0.227 0.258 0.348 0.212 0.197 0.333 0.242 0.197 0.197 0.258 0.288 0.288 0.258
0.303 0.242 0.258 0.242 0.212 0.333 0.273 0.318 0.288 0.303 0.273 0.212 0.182 0.288 0.303]
profilemostprobable(text, k, profile)



#GreedyMotifSearch(Dna, k, t)
#    BestMotifs ← motif matrix formed by first k-mers in each string from Dna
#    for each k-mer Motif in the first string from Dna
#        Motif1 ← Motif
#        for i = 2 to t
#            form Profile from motifs Motif1, …, Motifi - 1
#            Motifi ← Profile-most probable k-mer in the i-th string in Dna
#        Motifs ← (Motif1, …, Motift)
#        if Score(Motifs) < Score(BestMotifs)
#            BestMotifs ← Motifs
#    return BestMotifs


k=3 
t=5
dna=[
"GGCGTTCAGGCA",
"AAGAATCAGTCA",
"CAAGGAGTTCGC",
"CACGTCAATCAC",
"CAATAATATTCG",
]

function medianstring(motifs)
    ms = permutedims(hcat([[x for x in c] for c in motifs]...),(2,1))
    counts = vcat([sum([Int[x == y for x in ms] for y in ['A','C','G','T']][i],dims=1) for i=1:4]...)
    string([Dict(1=>'A',2=>'C',3=>'G',4=>'T')[x[1]] for x in findmax(counts,dims=1)[2]]...)
end

function score(motifs)
    stringform = [string(motifs[i,:]...) for i=1:size(motifs,1)]
    med = medianstring(stringform)
    sum(map(x->hamming(x,med),stringform))
end

function getprofile(motifs)
    stringform = [string(motifs[i,:]...) for i=1:size(motifs,1)]
    ms = permutedims(hcat([[x for x in c] for c in stringform]...),(2,1))
    vcat([sum([Int[x == y for x in ms] for y in ['A','C','G','T']][i],dims=1) for i=1:4]...)./size(ms,1)
end

function getprofilewithlaplacecorrection(motifs)
    stringform = [string(motifs[i,:]...) for i=1:size(motifs,1)]
    ms = permutedims(hcat([[x for x in c] for c in stringform]...),(2,1))
    vcat([1 .+ sum([Int[x == y for x in ms] for y in ['A','C','G','T']][i],dims=1) for i=1:4]...)./(size(ms,1)+4)
end

function greedymotifsearch(k, t, dna)
    bestmotifs = [d[1:k] for d in dna]
    bestscore = score(bestmotifs)
    for s=1:length(dna[1])-k+1
        motifs = String[]
        motif1=dna[1][s:s+k-1]
        push!(motifs,motif1)
        for i = 2:t
            profile = getprofile(motifs)
            motifi = profilemostprobable(dna[i], k, profile)
            push!(motifs, motifi)
        end
        newscore = score(motifs)
        if newscore < bestscore
            bestmotifs = motifs
            bestscore = newscore
        end
    end
    return bestmotifs
end


k=12
t=25
dna = [
"TTCGTCTTCGGCTGTATTGGAGTCAGACTGTGGAAGAATCCTAAGGCGGACATATAAACAAAGCCAATAAGAGTCCGTGACATGTCACGTCTATTAAGAGGCCTCGGGGTAAGATGCGCACGGTGGAAGTATCGTCAACAGCTGACGCCAAACGGG",
"CAAAAGAGTAGATTGATCTAATATTTAACTTATTTATATTAGGAGCCGCAGCTATGTGAAGTGTATCGAGACAGCGACTTCGACAGGCTGTCGGCATAACGTCCTCATGGTACCGAGTGCTTTTCCTGGTGAGGCGATGTAATAGTGATTTCACGG",
"GAATTAGTGAGGTCAACACCAAATGGACAACATGGAGCGGGCGGCTGCGGGTGCGACATGTCGGTCGAACAGGCCCAACTTACTAACAAGAGTAGCACGTGACCATCACGTATCGACTGTACATGTTGACATAACGTCTTCGACCGCTTCGCAAAA",
"AGCACCGTAACCTTCGTCTTCGCCCTTGTCCATGCTGACATGACAATCTGATCGGGCCTGGCCATTGATTCGACTCTGCAATAGCTTTGTTCGCGGGCCGTAAGTCCACGCTATCACGTCCACACTTGGCTCTCTGTTAATACGTGGTGACCACTT",
"TCTCCGACACGACCGGCAGCATGTTTAACGCTAACTTGAATGCTTCCCAACTCGTTCTGTTCCCTGCTGGTATCATATCGGTACATAGCACTACAATGCGCCTTCGACTGGAAAGGGGCGTGGACCCGATTCCTTCGTAGTATTGGCCAGGCCTGG",
"GGAGGCCTTCGTGTTTTCGCCCGCCCCGCGACAGCCAGGACGGTGTGCAGTCTCTCTCGTTCATACTTTACCACACATTCCTTTCATAGGGCGCGAACTCATCTGAAGCGCGTCTTCGACATTGCCGCGTCATTTCGGCATACACTCAAATAACGA",
"AACGTCTTCGACTATTGCTATCCATTAGGACGATCGAGTGCATCTATGCACCTAGGATAGGATCGCGTTCCCACCGTACTATCAATATGTTTGCCGCCTTCTTAAGCGCGCTCGTTGGTAGCGTTCTACTTCCAACTACGGTCTACAGCCAGCTCA",
"GCACCGGATATACTGCGGGCCAATTGGAATCTACGCGCTGTATCTATCTTCATTTGTCCATAAGTCTCACCCATTTAATGAACATTCGGCTTCGACCTCGAGTTTCTATTTATTGTACGGAGGGGCCCCACTCAAAAATACCCCTACGTTGGCCTT",
"CCACGCAACTTGTACGGCTTCGTCTCACGTAGACCAGATAGTAGCTTACCTGTGTAAGCCATACAGAATTCAAATTCCCGTCGTCAGCTGGTGCGGCGTAATTTCTACCATATGCGTTCTAGCAACTTAGCGCTAGACGAAGCCGTGATTTCGCGG",
"CTCAAACGCGATATTACCCCTCACCCAACGAGTTCTTACGTCTTCGTCACAGACCATCTAGGACGGATTCCGCCAAAGAGTCGTGTTTTAAAGGCATCGTCATCGCAGCTGAAGTTGTCCAGAGATGATATTGACAACTACATCCCTGACTCCTCG",
"AGTCGACGCCCGACCTGCAGGTACAGAGGCTAGCTCTTCGGCTTCGCCTACGATGCTTAAACGACCTCCCTCCCAGGAACGTTCGACCCTCATGAGAAACACCGCGTCATCCCTCAATAAGCACCTCCTTCAACGATTCGGTGCGATAGACATTCT",
"AATGCCGCCGATACTGAAGCGGCTCGTCTATGCTCCTAAGCAATTTGTATATATAGGCTTGTCTATTGGATACACGACTTCGTCAAGGACAACACGCGGTCAGAAAAGATCCGTCTTGTCACGCGGATTCCAGCGGTGTTTTAGCATCGTCTGAGC",
"AGCGCCTTCGGCGTTACAGGGGCCCGATTGCATCTTAGCCTCGGCATATCAGAAATGGGGTAATGATGAATGGCGTACTTAATTATGTAGGCGGTTCAAGTGTGGTCCTGCTCGTGCTGATGACGGAGGTGGGGCGGAAAATGGTGATTCGTCTCT",
"GTGTCGCGATTTCGTTTGAAGCACCATCGGCGCTTCGTTTGTGTCGCATCTGCCGGGAACATCGTCTTCGCCGACAGGATTAGATCATAGGGACCCGTCATCCACTCTTTGACTCCCACTCAATGCTCCACAGATTTCGTAGGTGTCTCAGCCTGC",
"TCCCTCAAGGTAAAATGATATGTTAACGCCTTCGGCGCACGTCTGCACTTCTATTATTAACAAGCCCATAGTGACGACAAATGAAACTAGTTTTAACATCACCCTTATCAAGTTAGTCCGCGTATTCCGGTACCCACAATGCTGTGGATATTCCCC",
"TGTGGGAGACGGGTCGCCTTCGTCAAGAAGCGGAACTTTTACTGAGCCGACGTGTTATCAAAGACCTTCAGATAGGTAAGAGAAGTGGCCTCACACTAACCTCATCAAAAGGAGACTTAAGCAATCCATTTACGTCCGAATTCGACATTGGTTTCT",
"ACGACACCCGCTGGCCGGTTTTCCATACCTCTCTTAAGTTACGATTCCGCGAGCCGCACCTGCGGCTTCGCCCGCCGCATACCAGGGAATTTAAGGATAAGCCGTGCGCCCACTGTCAGGTCGTTCCAGGAGAGCTTACAGCACTCCTCAGTATGC",
"ATGAGGCCAATTAACTATATATCTACTCGCAGGTCCGTCCAACGGCGCGGTCTACCGATTGCACTGAAAAGCTTGTTGAGTATGTTTCACCACGAGAAAGAGATCCCTATCCCGTGTGTTTCCGTCTTCGTCCGATCAACGTTCCCCCCACCTACG",
"GAATGATGATCGGAATACGTCCCGGAGCCGTTTCAAGCAATGTCCGTCGGGGTCTTAAAACGCTATTGACTCTTTTGTCCGGCCTTCGAGACCCGGTTCGGCTTCGCCCGCGCGTCTAGTAAAAGTAGTACCCCGTTCTACCTATCAGCCCTGTTA",
"TATCATTTTTTAGTCCTTACTAATGCTAGAAGAGGTCGATAACTGCACAGGAGCAGGCGGCAGGATCACCGTATTCTCCTTCCCGACCATGTTCCGCATCTGAGGAGGTCACATAGAACCCAAACCACGGTATGATGATAGTAAGTCGTCTTCGCC",
"GAGCCTAGTAACCCGCCTCCTAAACCGCTTCGACTAGGGTGCGTTGTTGGTATTACACCGGACGCCTTCGACCGACTATGACTACGCCTACCAGCTTTTTTCCTGACTACAGACTAAGCGCAACGTATATCAACTCCGGTGCCGCCGGATTCATGT",
"CAAAGTCGTTGTAGATGACAACACGAAAGACTGGGACACAGCTTGGCGAAGCCAAAGTATCAGTATCGTCGTGTAAAGGAGGCTGTCGCTATGATGATTTGGACACACGTCGGCTTCGGCGAGAAGTTTTCCAGTTCTCAATCCTAGGAAAAGCGA",
"TGTACGAGCCGTCCAACTTTAACGCTTCACCGGGGTTGTGACAGCTCAATGTGGAGAAATTCATACTCCGTTGGCGGCTTCGTCAGTTACTTGGCATACTTCCGAACTTTTAGACGGCTTCGCAACGATCGTGGAGCCCACAAAGAATCACAACGT",
"CCGCCGGTTCTCGGCTATCTAAGGGGCGCCTTCGCCTCAGTAGTAATTAAAGTGCACCACCGCAGGGCCTAAGAGAGGTAAGCTGTTCATGCCGGGTACTATCGGTGTAACTCCGCTATGGAATGAGTTCTAGCCTTGATCATAGATTGTACGTAT",
"GCTGACGTTCATATCTAGGGTATTTAAATCTTCACACTTTAAAACCTTACGCGACTAAGAGTTAGCATGCCGGATCTCCACGTTAAAGGGGCTTCCGATCTACATGTAGCATGACATTAGGCGATCAGCGTCCTTCGCCAGCATTGCGCCTTCGAC",
]
@time res = greedymotifsearch(k,t,dna)

function greedymotifsearch2(k, t, dna)
    bestmotifs = [d[1:k] for d in dna]
    bestscore = score(bestmotifs)
    for s=1:length(dna[1])-k+1
        motifs = String[]
        motif1=dna[1][s:s+k-1]
        push!(motifs,motif1)
        for i = 2:t
            profile = getprofilewithlaplacecorrection(motifs)
            motifi = profilemostprobable(dna[i], k, profile)
            push!(motifs, motifi)
        end
        newscore = score(motifs)
        if newscore < bestscore
            bestmotifs = motifs
            bestscore = newscore
        end
    end
    return bestmotifs
end

k=12
t=25
dna=[
"AGTGGGGTCCACGAGTATCAATTGAGGACCGTATATTCACAGGGCGTCGCAACAAAGACGACACATACTTTTACCCGTCATTCGAATTGCGTACACCTCGTGTGTGTCTCTGAGGTACGCGGAAAATTAATTTCAAGGGTCAGTATACGCTGATAA",
"TATGAACCGTCTAGTATCCGCAGGCTCCTTGAGCATGCCGGCTAACGTTACGCGATGTGGACATGCAACGTCCGCCCCCATCGATGAAACTACAAAGATTTCTAACATCAGCGAATAGAAGCGCTAAAGACGGAAGCTGCTAAGACCGTAACGGTA",
"CACGGGCGCACCGGTGCATGCTATGGCAGGACGTAGGGGTGGGCCCGCCTTTGACCTCTTGCATCAAAGACCTGCCGTCTAGCATTTTACAGGGTAGCGCTAGCGATAGAGGTCAGTTAATATGGCAGTTCCCCGCGATGAGGCTACCCGCCATGC",
"CTTCAAACACGGAGTGTTAGGTAAACCGATAGTGTGCCGTAAGTTCAGTCGGACCCGGTCCCAGGACCCGGCTACGCCGAGAAGTATCTCAGAACTTCTATCCGAGCATTGTTCAGTGTAGCTGGACACGCAGAGTATATCCGTGCATGAAAGACA",
"GTTACCGTACACTAACGATGGAACTGCCAGTACTCACTACGCCATTTGTCTCAGTGGGGTGAGACATCGTCCTGACTGCCGAACGCGAAAAAGACTGCTCTTTCGTCCCGACCTTAGCGCGAACATAAACACTTTTCTTTCTACCGGGGGTTCCTG",
"CCAAGGAGACAGACAAAACTTTGATTCAAGACGTACACTCTTGGCCGAGAACCAACACTATGGTTATTCTGTGCACAAAAGACAGGACCAGGCGGTATATGTCAGGGCTAATGTAGTACAGGGGACGGATGTCAAGGTGGCTCCGAGATAGAGTCT",
"TTCCTCGTTAATCTCTTACGGCTTGAAGTGCAGGCTAAATACAGAGTGGAGGCTCAATGCAGCGTATCGATTTTCAGCCGTGGTGCGGCAAAGACTGCGGCGTTGTCGGTTACAACGATACCATAGAGCTACCTGGAGGTCCCGGCCGAGGACCAC",
"GTTACCATGGAACGCATTTTTCCTTCTGGCAGCATTGCCTAATCTACCGCGCGAAAGACGTTTGGGTTCTGTTGAGCATCAGTATGCTAACACCACCCTTAGCTAGTTTCTTCCTTTTAGTACGTACCCTCCAATACTCTGTAATGCATCGCCCCA",
"AAAGGGTGATAATCGCGTCAGAGCAATCGCGTATCTTGAGACATAGGAACAGGGGCCTCCGAACTCCATGGTGATACGGCTGCTGCAACAAAGACTACTGAACGCTATGGTTAAAGCTGGCTTCACTGTCATCTGCTCGGTTGTGTCTGAAATGTC",
"GCCCAAAAGACAGCATTAGGGGCTATTTAGGAACGTAAGTCTACCTATTGTTAGGTCCCTTAGGCGCATTTAGGAGTGGAGTAGCATGGCCTGCTTGGCTGCAGGCTGATGGCTTGGGGAATGGATCACATCGCAGATGCCTAGGTATGGAATGTA",
"ATTGTACTTCAATCCTTCGGCCTTCGAACGTCGTCAGGGTAGTTGGTCGTGTTTACGACTATGCTAATGGCCATGCCCAGGTGTGGTATGGCAGCGAAGTCGTCGATGACGCTACTCGCCGCGTCAAAGACTAAAACTCTATTGCAACGAAGTTGC",
"ATTTGCAACATACTTTAGGACTGTGCGTCAAAGACGGCAAGGTTGAGATACTGTACAAAAAATTTTAGGGTAGGACGTATGGTCGGCCGTTTTCAATCTCTCATTAATTACTAAGTGCAAAAGCAAAATTTAAGTCAGGTGGAACAACGACACTAG",
"GCGACCAGCCTAGAAGGCCTGGCCAGGTGGCTAACTACAGTGAGGGTAGACTAGCCCATATAGTATGCGCCTCATCGACACAACCGACCATGCCGCTTACACCAGATAGCACGAAAGACGCGTCTCTTTGCGCAATTGAAAAAGCATCCGACGTGT",
"CTTCAGGGTGCAGGTGAGAAAATACGTGGGAAGCTCAGAGCCACTGCCCCCTCGCCGTGGGCAGGAAGCTCCTAACCGCAGGGAGCACCAAAGACCAGGTCCGCGGGAGGCGCCACTCAGTGTATCTGCGTATTGTTTAGGAATATAGGAATGTCC",
"ATGCAAGCGGCAGTCCTGTCCGTTCCTTCAAATATTGATCTTTCGAACGGGTTAAGAAGACGAGGAAGCCATGATCCATGTTAGCCGGGGCTCGCAGGGTGCACCCCTGTGTCCATGCCTGCCCGAAAGACAGGGCGTAGTCAGGCTTTACGCGAA",
"ATACTATCCCCACACGGATGCTGCGCGGTAAAGACTGCACAGCTGCGCATCTGTCATCATAGTGCTAATAGAGGCAATACGATGTGATGGCGTCGGTCGAATTGAAACTCATTGCAATTCGCGCCATACGTGGATCACTCGGCAACAAGGAGCAGT",
"CCGATCTCCGCCACTGGCGTGGACGGTTGCAATAAAGTAAAGTTCCGGTCGGCTGTACGTTAGAAAGATAAAGGGTTCGTGCTTAACGATATGCTCAACCGGGTCTACCTACGTGGTGGCGCGTCAAAGACGACTTCTGCCGGAGAGAAGCAAACG",
"ACGAGTGTGTGGAATGAAGACATGTCCTAGATGACTGACCCTGCATAATGCACGGAAAGACACGAATCCACCTCGGGATGTCGACAAAAATACTTTGCGACAAAGACAGACACACGAAAAATAACTTCTACCGGGACGTACCCTCGCAGGGCACAG",
"CAAGAATCGTAGTTCTCCTCTAACGGCCTAGAAATTTCGCCATTTGTATTTATCAGGTAAAGTAAGCAGTACATTTATTACTCGCAATCAACGCAAAAGGACCCCTCCTCCCCGCGTAAGCATTGGTGCTCTGCGTCAAAGACACTCTGCATCGAC",
"GAACTCGCTACGAAATGGTCCCTGTATGACCCTGGTTTTCATTGTTGTCGGTCGCTGGATGCTGGAAAGACGTTGACGGTTAACAGCTAGGTAAAACAAATCGTGGGCGGCGCTGCATCCTTAATCGATTTTTTTTGGCAAGGCTATGGCGTCGGG",
"TCGCCCAAACGTATATCGGCAGCTTTCTGCCCCGCAGTACAAGTCACGAAAAAAATCGCGAGGTAGCGGCACCGCCTTCAAGTCGATGTCATTGTCAATCGAATCGTGGGCATTATGTCAGCGTCAAAGACCATCACGTACAAGACCTGAACGACC",
"TTCGCCGCCTTGTTTCGAAGTTGGCCGTCGGGGTGAGACAAGACATTGGCACCTAGCGAACCCATGACCGTCAGTCCCCTGATAGCCGTAAAGACTCCGGTGACAAAAAACGAGGGATGCCCACTGTGCTGAATCATCCGAGAGCACGTCTCGGTG",
"CCAGCGAGTATCTATCCCAACTGATTATAAGGGTGAGCGTAAAAGACGACTACCCGACTTCTTAACTGGTTTCAAGTATATCTTGTCCTCCGCACGCATTCACCACACGTCCGTCTAAACTTCCGACAAATGTCTTAATGTCCCATTAACCTGGCG",
"CACGTTAAGACTTAAGCCTGTGTTGTAAACTGATATGCAGAAAAGACCACGACGATATCCACGGTACGATGTTCTGCGGTGCGAACTTCCGCTTAAGGTTAGATCCCGTGCTAGTTGTGTTACTTGAAATAATATCCCTCCCGTAGTTAGCAAGAA",
"TCCCAGGTAAGTTCGACATACGGTACTTCAGAGATCCGCTGGTGGCGTTGACCAGAAGGTCGTGTCTAAAAAGGGGTGCGCAGCGCAATCGTCACACCAGGCATGATTTAAACGAAGTGGGCTAGAAAGACATGGGTGGGACTTATTAAAGTCATT"
]
greedymotifsearch2(k, t, dna)

function randomizedmotifsearch1(dna, k, t)
	motifs = [begin; i=rand(1:length(s)-k+1); s[i:i+k-1]; end for s in dna]
	bestmotifs = motifs
	while true
		profile = getprofilewithlaplacecorrection(motifs)
		motifs = [profilemostprobable(s, k, profile) for s in dna]
		if score(motifs) < score(bestmotifs)
			bestmotifs = motifs
		else
			return bestmotifs
		end
	end
end

function randomizedmotifsearch(dna, k, t, n=1000)
	motifs = randomizedmotifsearch1(dna, k, t)
	bestmotifs = motifs
	for i=1:n-1
		if i%1000==0
			println(i,bestmotifs)
		end
		motifs = randomizedmotifsearch1(dna, k, t)
		if score(motifs) < score(bestmotifs)
			bestmotifs = motifs
		end
	end
	return bestmotifs
end


k=8
t=5
dna=[
"CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]

randomizedmotifsearch(dna, k, t)

k=15
t=20
dna=[
"ATAGCTCTATCGTACGTCTATACATCAAGCTACCGAGATAGGGTGCTGGCGTGCCAGCCGGCTGTCCTGGAGACATATTGGGCCGCGAAACGCCACAGACTACACTGGGCCACCAACGGCTTACCGTGAGCACGCATAAACCTGGAAGGGATAGCTCTATCGTAC",
"GTCTATACATCAAGCTACCGAGATAGGGTGCTGGCGTGCCAGCCGGCTGTCCTGGAGACATATTGGGCCGCGAAACGCCACAGACTACACTGGGCCACCAACGGCTTACCGTGAGCACGCATAAACCTACGTAGCAGAGCAGTGGAAGGGATAGCTCTATCGTAC",
"CAATGACTGGTAGATCTTACTACACCTGACACTGTACTTACTGCGCAACACAGCTGGACTCACATTTATGATTGTAGAGAGGTCATAGTAATTGGCCAACGTCCCACACTCACTGTTGTTAAGCTAATTCTCCGACCGAGCGGCTTAGCGCAGTCTATCGTAATT",
"AAGAATTGAAACGTAACAATTTGTCAAGCTTTTCTCCTAGATTACTAGTTCTTCGGATCAGCTTGTCACTTCGGGCCTATGGTATGCGCAGAGCAGGTCCATTATCGCCGCGACGCCGGTACGACCTAAAGACTTTTTATGTATCGGGTATTACCCACGCCATCT",
"CCGACAGACCTGGCCCAAATCGAATCTCGTGGACGAAGTTGGGTATCCATTCGTCCCATTCTTTTGTTAGATGCTAGCTACATACGAGCACACGTTTCCGTGTCTTATCACTTGCTTGGTAAAATCGGGGTACTGCGCAGACGCGGCGGATTGAATCCCCATACT",
"ACAGTGGTTTCTATTGCTCCTCTTTTTATACTGCAGTGAGCAGCGTTCATAGCTGCGATTTCGCCTTGGCCGGGCCCTATACACATGTAGTGCCCAGGCATTGGTCCGGCTTCAGATCATCGCTTCATAGAAATGTATCGCATTAATAGGTATGTAGCCCAAGCC",
"GGTGGCTCTGTTTTCGGTGTGGATTGCGCCGCCTCCTACTGCGGGTAGCAGGGCCCCCACTGTAACCATTTCGCTTGGTGGAGATCTTAATCAATTGAAGGACGCTGCCTCCGATAGTTATCCGTAGCGCTTGTGACGAAGGGCACTCATGCGATGCGATGGCAA",
"GGTAACTTTTGTCGTGGGATGTAGGGTGGCTTACACCGGTACAGAATTACGCGCAGAGCAGCTCGAGCTTAGGTATCTTAGGGAACACCGGAAAACTTAAGTCAGACAAGTGGGTCACCTTGGCGACCTTAGGTGGATGATGGACAAGACTCCCCTTGCCATTGA",
"ATTTTCAATTCCAACTCTGACTCCATGTATAGAATTCATAACCCAGAGCCTATGTACTATTCACTGTGTGCCTTGTAATCAATGGCTCGAATGACTGCACTGCTGCCTCTCGTGGAGTCAGGTGAATACCACCTTACTGCTGGGAGCAGGGTCGGGGCCATGGAA",
"AAAAATGCCTGCGATATGTGTCATTCTCTTAATCGACACAACTATACCGCCAAAATGCCTTCGAGCCAAGTCTCTCTTCCGTACTGCGCAGAGGCTCCATTAACGAGGGCGTGACGCGTACCTAGCCAACGGGGCGGAGCCCCGAGTCCAAGGGACCGCTGCTAC",
"TAAGACCTCAGGTCCCGCCGACTACAGAAGAATTATGAAGTAACATGCGGAAGGCGGCTGTAAGATGCGGTAGCACGCAGAGCAGAGTATGAAGCCAAGAACCACAACGTCTGCATTATTCTCATGCGCACCTGATCAGCGCCTGCCACTCTCTCGGTACCGGAG",
"TCATCCCCGCTTAGCTCACACTGTTGCCAATTTTCGTAAATGGATATAGATTACTGAGGCTGTTGGCTTTTAGAAGAATTATCACATTGAATATGGAAATATGTACTGCGCAGGAAAGACGGTATAGACGACTGTGGAGACGCGCCAGTGAAAGCGCGCCCTGGC",
"GTGGAGCATCCCGAATCTCCAACGTGCGGACTCACACCCATCGGTAATTCACTGGTACATCTGACGCAGAAGGCACATGCGTAAAGAGGTGCTGAGGGATGAAGATGCCGTGTGTACTGAATTGACCATCTGTGTGATAGGTACTGTAAAGAGCAGACAGGAATG",
"GGTGACTTAGGCAGGAGACAAAGTGTAACTTCGGATACATGGCAGAGCAGCCCTACTAATCGAAGGCCACAGCGATATAGGGTCCGATCTCTCGACACAGTCGTGGCAGTTGTTTGGAAATTATAGCTCTCTACAAGACGAGCGGCACGGCCCCCTCCGCGACCG",
"CTTATACGTCTGAGCCAAGATATGATTTTACTTGCACGATGCCTGGATGGGAAGCTGTCTAACTGCTAAACTCGGAGAGCCCCGTTGGAATCTAAAGAGCATCTCCCCCGACACCTACTTGATACACCACATCTGCGCAGAGCACTTCATAATCCTGTATCCCGA",
"GCCGAGCTAATGACTCACTATAGTGGTTGCGTATGGTATCCAATTGACATTGGTAGACTGCGCAGAGCGTGGTTCCGGTGCCCCCGAGTGTGAACGACTTATGCTGCGAGACGTATAGCGGCTTCTATCTGTATACTGCATGGGGGAGATACACGACCCGCATGC",
"AGAAGCCCATTTTAGATGACCGGGTCGTCTTGGATCAACTACTAGGGATCACTTGGATCGGTCGTTCACTGACGTTGTCCCCGTCCGGAAAGCGTATTACTGCGCCAGGCAGACGCTGGTTAGTGAGATCACCTCTCGGCAGCTGCTTAGCAGAGAGTGAGACCG",
"AGCGCTTTTCGGGTGAGGCAAACATATTGCGCGTCGCCGATGTTCGTTTCACTTTCATCTACTATTGCTCATCACGATCCCAAGATACTTGTTAGCCTCCTAGATAGCATGGAGCGCGTCATTTCCGATGTCCTCGGTACTAGCCAGAGCAGTTCCGCGTTCCCC",
"CAGACTAGCAGGATCCTCCCCGGCGGAAGGTACTTTTCAGAGCAGTGGAGGTACATTCTTAGTGCCCAAAGCTCTGTATATCCGCATGAGAATGGTTGAGTCAAGGAGGCTCAATAGCTACAGCTGTCAAGACTACTATGTAGACCTTGGAACCGCACCTTGTAC",
"AGAACTTACTGATAAGAGCAGGAAGGTCCATGATAGGTGGGGCTATCCAACCTTTCTCGATGAGTTGCTAATGATACTTATATAGGACCACCAATTCTTGCACGTCATCCGGGATTGAGTATAGTGAGCTTACCAGCCCGAGGTTATCATAGCGAGGCCGGCTAG"]

randomizedmotifsearch(dna, k, t)


#Code Challenge: Implement GibbsSampler.
#Input: Integers k, t, and N, followed by a collection of strings Dna.
#Output: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts. Remember to use pseudocounts!
#Note: The next lesson features a detailed example of GibbsSampler, so you may wish to return to this problem later.

#GibbsSampler(Dna, k, t, N)
#        randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
#        BestMotifs ← Motifs
#        for j ← 1 to N
#            i ← Random(t)
#            Profile ← profile matrix constructed from all strings in Motifs except for Motifi
#            Motifi ← Profile-randomly generated k-mer in the i-th sequence
#            if Score(Motifs) < Score(BestMotifs)
#                BestMotifs ← Motifs
#        return BestMotifs

function gibbssampler(k,t,N,dna)
	motifs = [begin; i=rand(1:length(s)-k+1); s[i:i+k-1]; end for s in dna]
	bestmotifs = motifs
    for j = 1:N
        i=rand(1:t)
        getprofilewithlaplacecorrection(vcat(motifs[1:i-1], motifs[i+1:end]))
        ...
    end	
end

function randomizedmotifsearch1(dna, k, t)
	motifs = [begin; i=rand(1:length(s)-k+1); s[i:i+k-1]; end for s in dna]
	bestmotifs = motifs
	while true
		profile = getprofilewithlaplacecorrection(motifs)
		motifs = [profilemostprobable(s, k, profile) for s in dna]
		if score(motifs) < score(bestmotifs)
			bestmotifs = motifs
		else
			return bestmotifs
		end
	end
end

function randomizedmotifsearch(dna, k, t, n=1000)
	motifs = randomizedmotifsearch1(dna, k, t)
	bestmotifs = motifs
	for i=1:n-1
		if i%1000==0
			println(i,bestmotifs)
		end
		motifs = randomizedmotifsearch1(dna, k, t)
		if score(motifs) < score(bestmotifs)
			bestmotifs = motifs
		end
	end
	return bestmotifs
end

k=8
t=5
N=100
dna=[
"CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

#TCTCGGGG
#CCAAGGTG
#TACAGGCG
#TTCAGGTG
#TCCACGTG

gibbssampler(k,t,N,dna)

end # module
