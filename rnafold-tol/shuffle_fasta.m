% Read a (possibly large) fasta file; shuffly each sequence using shuffleCDS(); write the output to another fasta file.
% TODO: geneticCode should be specified for each sequence, since different sequences from the same organism are translated using different codes.
function []=shuffle_fasta(filename,geneticCode)
  % Read sequences from fasta file; write shuffled sequences to output file.
  biofile = BioIndexedFile('FASTA', filename);
  numSequences = biofile.NumEntries;

  outfile = [filename '.shuffled.v2.fna'];

  for i=1:numSequences
    entry = biofile.read([i]);

    proteinId = entry.Header;
    seqBefore = entry.Sequence;
    seqAfter = shuffleCDS(seqBefore, geneticCode);

    fastawrite(outfile, proteinId, seqAfter);
  end

end
