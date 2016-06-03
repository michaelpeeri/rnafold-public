% Shuffle a CDS sequence while maintaing translation and codon frequencies.
%
% This function simply randomly permutes all synonymous codons within the CDS.
% (maintaining 1. The translation and 2. The codon frequencies (and related
% compositional characteristics such as GC content)).
% Using this type of shuffling simulates a situation in which those are the
% only characteristics under selection.
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
  origCdsCodons = cdsCodons;

  lenInCodons = size(cdsCodons,1);

  % Perform shuffling for occurences of a given AA
  function ret = shuffleAA(aa, aaCodons)
    numReplacementsDone = 0;
    ret = struct('aa',aa, 'replacements', 0);

    if(strcmp(aa, 'Stops') || strcmp(aa, 'Starts') || strcmp(aa, 'Name'))
      return;
    end

    % Skip AAs not encoded by multiple codons
    if( size(aaCodons) < 2 )
      return;
    end

    occurrences = find(translation==aa);     % occurrences is a boolean vector indicating the positions of aa

    % Skip AAs which occur only once
    if( length(occurrences) < 2)
      return;
    end

    %disp(['------------------[ ' aa ' ]--------------------']);
    %disp(aaCodons);


    selectedCodons = cdsCodons(occurrences,:);

    cdsCodons(occurrences,:) = selectedCodons(randperm(size(selectedCodons,1)),:);


  end
  % Perform shuffling individually for each AA
  x = cellfun(@(aa) shuffleAA(aa, getfield(rcode,aa)), fieldnames(rcode));
  %disp(struct2table(x));

  numEqualCodons = sum(all(cdsCodons == origCdsCodons, 2));
  codonReplacementFraction = (lenInCodons - numEqualCodons) / lenInCodons;
  assert(codonReplacementFraction >= 0.0 && codonReplacementFraction <= 1.0);
  if( codonReplacementFraction < 0.3 )
    disp(sprintf('Warning: only %d replacements made (%2.1f%%)', (lenInCodons - numEqualCodons), codonReplacementFraction*100.0));
  end

  newcds = reshape(cdsCodons', 1, []);
  if(~strcmp(translation, nt2aa(newcds, 'ACGTOnly', false)))
    disp('Warning: Translation might have changed!')
    disp(translation);
    disp(nt2aa(newcds, 'ACGTOnly', false));
    % Ignore this (for the moment!)
  end
end
