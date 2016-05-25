
function [] = process_genome(fromEntry,toEntry)
  % Configuration
  %filesToProcess = struct('path',{},'geneticcode',{}, 'id',{});
  %filesToProcess(1).id = 'Chlamydomonas reinhardtii';
  %filesToProcess(1).path = 'GCF_000002595.1_v3.0_rna.fna ';
  %filesToProcess(1).geneticcode = 1;
  %filesToProcess(2).id = 'Phaeodactylum tricornutum CCAP 1055/1';
  %filesToProcess(2).path = 'GCF_000150955.2_ASM15095v2_rna.fna';
  %filesToProcess(2).geneticcode = 1;

  % read genome sequences into indexed object
  % TODO - receive these as parameters...
  %genomeFile = '~/rnafold/data/GCF_000002595.1_v3.0_genomic.fna';
  %gffFile = '~/rnafold/data/GCF_000002595.1_v3.0_genomic.gff';
  %species = 'Chlamydomonas reinhardtii';
  %geneticCode = 1;
  genomeFile = '~/rnafold/data/GCF_000150955.2_ASM15095v2_genomic.fna';
  gffFile = '~/rnafold/data/GCF_000150955.2_ASM15095v2_genomic.gff';
  species = 'Phaeodactylum tricornutum CCAP 1055/1';
  genomeSeqs = BioIndexedFile('fasta', genomeFile);
  geneticCode = 1;  % ref: http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=556484


  results = fopen(sprintf('~/rnafold/results_%d_%d.json',fromEntry,toEntry), 'w');
  fprintf(results, '{\n "Species":"%s",\n "geneticCode":%d,\n "genomeFile":"%s",\n "gffFile":"%s",\n "script":"%s",\n "time":"%s",\n "CDSlist":[\n', species, geneticCode, genomeFile, gffFile, mfilename('fullpath'), datestr(now()));

  % return substring from genoming sequence
  function seq = getGenomeSubsequence(id, from, to, strand)
    entry = fastaread(getEntryByKey(genomeSeqs, {id}));
    seq = entry.Sequence(from:to);
    if( strand=='-')
      seq = seqrcomplement(seq);
    end
  end
  % test
  %disp(getGenomeSubsequence('NW_001843466.1', 10, 20));

  % helper - format value according to type
  function s = formatCell(x)
    if( isnumeric(x))
      s = num2str(x);
    elseif( islogical(x))
      if(x)
        s = 'true';
      else
        s = 'false';
      end
    else
      s = ['"' x '"'];
    end
  end

  % Parse the genome annotation (contained in a GFF-format file)
  gff = GFFAnnotation(gffFile);
  % Filter the features by type 'CDS' only (ignore the rest of the features)
  gffCDS = getSubset(gff, 'Feature', 'CDS');
  %a = getData(gffCDS, [1:100]);

  % Helper function - convert the CDS 'Attributes' value to a matlab struct
  function atts = parseCDSattributs(attsStr)
    keysAndValues = strsplit(attsStr, {';','='});
    atts = struct(keysAndValues{:});
  end

  % Process all CDS entries
  count = 0;
  prevCDS=struct();
  prevatts=struct('ID','');
  sequence = '';
  numExons = 1;
  firstEntry = true;

  %
  % Sample CDS entry:
  %
  % Reference: 'NW_001844015.1'
  % Start: 293
  % Stop: 343
  % Feature: 'CDS'
  % Source: 'RefSeq'
  % Score: '.'
  % Strand: '-'
  % Frame: '1'
  %
  % Sample attributes:
  %
  % ID: 'cds228'
  % Parent: 'rna228'
  % Dbxref: 'InterPro:IPR000626,Genbank:XP_001702404.1,GeneID:5727945'
  % Name: 'XP_001702404.1'
  % Note: 'post-translationally modifies proteins to tag them for proteasome-mediated degradation%2C internaliza…'
  % end_range: '5927,.'
  % gbkey: 'CDS'
  % gene: 'UBQ3'
  % partial: 'true'
  % product: 'bi-ubiquitin'
  % protein_id: 'XP_001702404.1'
  %

  function processCDSEntry( CDS )
    %disp(CDS);
    atts = parseCDSattributs(CDS.Attributes);

    % Get the sequence for this CDS
    if(~isfield(CDS, 'dummy'))
      newseq = getGenomeSubsequence(CDS.Reference, CDS.Start, CDS.Stop, CDS.Strand);
    else
      newseq = '';
    end

    function emitCompletedCDS(CDS, atts, sequence)
      count = count + 1; % Count the number of CDS-containing genes found
      disp('-------------------');
      disp(atts.ID);
      disp(CDS);
      disp(atts);
      disp(sequence);
      disp(length(sequence));
      disp(numExons);

      if(~firstEntry)
        fprintf(results, '  ,');
      else
        fprintf(results, '  ');
      end
      attsJson = strjoin(cellfun(@(k,v) sprintf('"%s":%s',k,formatCell(v)), fieldnames(atts), struct2cell(atts), 'UniformOutput', false), ', ');
      fprintf(results, '{\n  "id":"%s",\n  "CDS":"%s",\n  "numExons":%d,\n  "length_nt":%d,\n  "attributes":{%s},\n  "reference":"%s",\n  "strand":"%s"\n  }\n', atts.ID, sequence, numExons, length(sequence), attsJson, CDS.Reference, CDS.Strand);


      if(isfield(atts, 'transl_table'))
        currGeneticCode = geneticcode(str2num(atts.transl_table));
        disp(['Overriding genetic code: ' atts.transl_table]);
      else
        currGeneticCode = geneticcode(geneticCode);
      end

      translation = nt2aa(sequence, 'GeneticCode', currGeneticCode, 'ACGTOnly', false);
      numStopCodons = sum(translation=='*');
      disp(translation);


      % Do sanity tests
      assert(length(sequence)>3);
      %assert(mod(length(sequence),3)==0);  % TODO - This doesn't work for partial CDSs - why?
      isPartialCDS = isfield(atts, 'partial') && strcmpi(atts.partial, 'true');
      if(~isPartialCDS)
        assert(mod(length(sequence),3)==0);
        if(~isfield(atts, 'transl_except'))
          assert(any(strcmpi(currGeneticCode.Starts, sequence(1:3)))) % Is the first a codon a start-codon?
          %assert(strcmpi(sequence(1:3), 'atg'));
          assert(numStopCodons==1); % A complete CDS must contain a stop codon (and it must be the last codon - see below)
        end
      end
      if(~isfield(atts, 'transl_except'))
        assert(numStopCodons<=1); % There is at most one stop codon

        if( numStopCodons > 0 )
          assert(translation(length(translation))=='*'); % if a stop codon exists, it must be the last codon
        end
      end
    end

    if( ~strcmp(atts.ID, prevatts.ID) && length(prevatts.ID)>0) % is this CDS not a continuations of the previous gene?
      % Emit the previous sequence
      emitCompletedCDS(prevCDS, prevatts, sequence);

      % Setup the next sequence
      sequence = ''; % reset the sequence buffer
      numExons = 1;
      firstEntry = false;
    else
      % Continuation
      numExons = numExons+1;
    end
    sequence = [sequence newseq];
    prevatts = atts;
    prevCDS = CDS;
  end
  arrayfun(@processCDSEntry, getData(gffCDS, [fromEntry:toEntry]));
  % Append a dummy CDS to emit the last actual CDS
  processCDSEntry( struct( 'Attributes', 'ID=dummy__', 'dummy', 'true'));

  disp(sprintf('Found %d entries containing CDSs', count));

  fprintf(results, ' ]\n}\n');
  fclose(results);




  %[Header, Sequence] = fastaread('~/rnafold/data/Scer_coding_regions_from_biomart.fasta');
  %seqs = fastaread('~/rnafold/data/GCF_000002595.1_v3.0_rna.fna');

  %e1 = arrayfun(@(s) getEnergy(s.Sequence), seqs(1:2));
  %disp(e1);

  %e0 = arrayfun(@(s) getRandomEnergy(s.Sequence), seqs(1:20));
  %disp(energies);
end

function codons=splitIntoCodons(seq)
  disp(length(seq));
  assert(mod(length(seq),3)==0);
  codons = cellstr(reshape(seq,3,[])');
end


function energy = getEnergy(seq)
  %[ignore, energy] = rnafold(seq);
  [status,cmdout] = system(sprintf('python2 ~/rnafold/rnafold_seq.py %s', seq));
  assert(status==0);
  energy = str2double(cmdout);
  disp(energy);
end

function energy = getRandomEnergy(seq)
  codons = splitIntoCodons(seq);
  disp(codons);
  swapCodons(codons, 3, 5);
  disp(codons);
end

function swapCodons(seq, idx1, idx2)
   tmp = seq(idx1);
   seq(idx1) = seq(idx2);
   seq(idx2) = tmp;
end
