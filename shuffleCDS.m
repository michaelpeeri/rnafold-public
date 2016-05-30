% Shuffle a CDS sequence while maintaing translation and codon frequencies.
% This function maximizes the number of replacements, by avoiding cases in which a codon is replaced by an identical codon.
%
% Explanation:
% Each AA is handled in sequence, with replacements done between the different codons used to encode it.
% The problem can be formalized as a graph, with each codon instance (encoding the current AA) represented as a node and each potential
% replacement as an edge. Edges exist therefore between a node and all nodes representing instance of other codons. This kind of graph
% is called a Complete multi-partite graph (multi-partite because it contains several classes of nodes with no edges within the same
% group; and complete because all edges to other groups exist). In this formalization, we are looking for (one of) the maximum-cardinality
% matching.
% The algorithm here is similar to (but not the same) as the one developed in this reference:
% David Sitton - Maximum Matchings in Complete Multipartite Graphs
% http://math.furman.edu/~mwoodard/fuejum/content/1996/1996paper2.pdf
% In my version, all possible replacements are done between the two most frequent codons (for the current AA), and this is repeated as
% long as possible replacements remain. The intuition behind choosing the most frequent codons for the first replacement (this is a
% greedy heuristic) is that these are the most constrained groups (because there are no edges inside a group, the largest group is the
% hardest to find matches for. A codon from a small group can be replaced with most other codons).
%
function newcds=shuffleCDS(cds, geneticCode)
  % TODO: use each sequence's genetic code
  % TODO: special handling for stop codons ?
  assert(mod(length(cds),3)==0);
  code = geneticcode(geneticCode);
  rcode = revgeneticcode(geneticCode);
  translation = nt2aa(cds, 'GeneticCode', geneticCode, 'ACGTOnly', false);
  cdsCodons = reshape(cds,3,[])'; % convert the sequence into a Nx3 matrix
  assert((size(cdsCodons,2)==3) && (size(cdsCodons,1)==length(translation)));


  % Perform shuffling for occurences of a given AA
  function ret = shuffleAA(aa, aaCodons)
    numReplacementsDone = 0;
    ret = struct('aa',aa, 'replacements', 0);

    if(strcmp(aa, 'Stops') || strcmp(aa, 'Starts') || strcmp(aa, 'Name'))
      return;
    end

    %disp(['------------------[ ' aa ' ]--------------------']);
    %disp(aaCodons);

    % Skip AAs not encoded by multiple codons
    if( size(aaCodons) < 2 )
      return;
    end

    occurrences = (translation==aa);
    % Skip AAs which occur only once
    if( sum(occurrences) < 2)
      return;
    end

    % Find all occurrences of each codon
    lenInCodons = size(cdsCodons,1);
    codonOccurrences = cell2mat(cellfun(@(x) all(cdsCodons==repmat(x, lenInCodons,1),2), aaCodons, 'UniformOutput', false));
    codonNumOccurrences = sum(codonOccurrences,1);

    % Helper function to find 'count' random positions in which a certain codon occurs
    function indices = getRandomCodonPositions(codonIdx, count)
      %disp(sprintf('-- getRandomCodonPositions(%d,%d)', codonIdx, count));
      allIndicesForCodon = find(codonOccurrences(:,codonIdx));
      %disp(allIndicesForCodon);
      randomized = allIndicesForCodon(randperm(size(allIndicesForCodon,1)));
      %disp(randomized);
      indices = randomized(1:count);
    end


    % Main loop - perform codon shuffles (for codons encoding for the current AA) as long as possible
    while(true)
      %disp('Sorting...');
      % Get the number of occurrences for each codon encoding for the current AA.
      % B = Counts; I = Indices
      [B I] = sort(codonNumOccurrences, 'descend');
      %disp(B);
      %disp(I);

      Disp(codonNumOccurrences);
      topCounts = B(1:2);
      topCodons = I(1:2);
      %disp(topCodons);
      %disp(topCounts);

      if( topCounts(2)==0 )
        disp('No more replacements are possible.');
        break
      end

      replacementCount = topCounts(2);
      codonA = topCodons(1);
      codonB = topCodons(2);
      %disp(sprintf('%d replacements possible between %d <-> %d', replacementCount, codonA, codonB));

      % find and index for codon A
      indicesForCodonA = getRandomCodonPositions(codonA, replacementCount);
      indicesForCodonB = getRandomCodonPositions(codonB, replacementCount);
      %disp('will replace using indices: ');
      %disp([[indicesForCodonA indicesForCodonB]]);

      ntOfCodonA = aaCodons{codonA};
      ntOfCodonB = aaCodons{codonB};
      whos ntOfCodonA;
      %disp([sprintf('Replacing %d occurrences of ', replacementCount)  ntOfCodonA '<->' ntOfCodonB]);

      %disp(size(codonOccurrences()));

      for i=1:replacementCount
        % Exchange the codons

        assert(strcmp(cdsCodons(indicesForCodonA(i),:), ntOfCodonA));
        assert(strcmp(cdsCodons(indicesForCodonB(i),:), ntOfCodonB));

        cdsCodons(indicesForCodonA(i),:) = ntOfCodonB;
        cdsCodons(indicesForCodonB(i),:) = ntOfCodonA;
        % Update the codon occurrences table
        assert(codonOccurrences(indicesForCodonA(i),codonA)==1);
        assert(codonOccurrences(indicesForCodonB(i),codonB)==1);
        codonOccurrences(indicesForCodonA(i),codonA) = 0;
        codonOccurrences(indicesForCodonB(i),codonB) = 0;
        assert(codonOccurrences(indicesForCodonA(i),codonA)==0);
        assert(codonOccurrences(indicesForCodonB(i),codonB)==0);
      end

      % Remove the occurrences just performed from the occurence counts vector
      codonNumOccurrences(topCodons(1)) = codonNumOccurrences(topCodons(1)) - topCounts(2);
      codonNumOccurrences(topCodons(2)) = codonNumOccurrences(topCodons(2)) - topCounts(2);

      numReplacementsDone = numReplacementsDone + replacementCount;

    end

    ret.replacements = numReplacementsDone;

  end
  % Perform shuffling individually for each AA
  x = cellfun(@(aa) shuffleAA(aa, getfield(rcode,aa)), fieldnames(rcode));
  %disp(struct2table(x));

  newcds = reshape(cdsCodons', 1, []);
  if(~strcmp(translation, nt2aa(newcds, 'ACGTOnly', false)))
    disp('Warning: Translation might have changed!')
    % Ignore this (for the moment!)
  end
end
