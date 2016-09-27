
function [] = extract_cds_from_gff(fromEntry,toEntry)
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
  genomeFile = '~/rnafold/data/JGI/PhytozomeV11/Creinhardtii/assembly//Creinhardtii_281_v5.0.fa';
  gffFile =    '~/rnafold/data/JGI/PhytozomeV11/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene.gff3';
  species = 'Chlamydomonas reinhardtii';
  %geneticCode = 1;
  %genomeFile = '~/rnafold/data/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa';
  %gffFile = '~/rnafold/data/GCF_000150955.2_ASM15095v2_genomic.gff';
  %gffFile = '~/rnafold/data/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff';
  %species = 'Saccharomyces cerevisiae'
  genomeSeqs = BioIndexedFile('fasta', genomeFile);
  geneticCode = 1;  % ref: http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=556484

  % Results fasta file
  resultsFasta = sprintf('~/rnafold/results_%d_%d.cds_only.fna', fromEntry, toEntry);
  if( size(dir(resultsFasta),1) == 1)
    disp(sprintf('Warning: overwriting fasta results file %s', resultsFasta));
    delete(resultsFasta);
  end

  % Results json file
  results = fopen(sprintf('~/rnafold/results_%d_%d.json',fromEntry,toEntry), 'w');
  fprintf(results, '{\n "Species":"%s",\n "geneticCode":%d,\n "genomeFile":"%s",\n "gffFile":"%s",\n "script":"%s",\n "time":"%s",\n "CDSlist":[\n', species, geneticCode, genomeFile, gffFile, mfilename('fullpath'), datestr(now()));

  % return substring from genoming sequence
  function seq = getGenomeSubsequence(id, from, to, strand)
    %disp(id);
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

    if(~isfield(atts, 'ID'))
      if(isfield(atts, 'Name'))
        atts.ID = atts.Name;
      end
    end
    % Remove trailing .N (e.g., ID.1, ID.2, ID.3), used in CDS identifiers in recent JGI GFF3s.
    pos = regexp(atts.ID, '[.]\d+$');
    if( pos )
      atts.ID = atts.ID(1:pos-1);
    end
    assert(isfield(CDS, 'Strand') || isfield(CDS, 'dummy'));
    assert(isfield(CDS, 'Start')  || isfield(CDS, 'dummy'));
    assert(isfield(CDS, 'Stop')   || isfield(CDS, 'dummy'));

    %disp(CDS);

    % Get the sequence for this CDS
    if( ~isfield(CDS, 'dummy'))
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
      %disp(sequence);
      disp(length(sequence));
      disp(numExons);

      isPartialCDS = isfield(atts, 'partial') && strcmpi(atts.partial, 'true');

      % Skip CDS that are not marked as 'Verified'
      % Todo - make this configurable, so it works with organisms other than S. cer.
      %
      %if( isfield(atts, 'orf_classification') && strcmpi(atts.orf_classification , 'Verified'))
      %if( ~isfield(atts, 'orf_classification') )
      %  disp('Skipping unclassified ORF');
      %  return;
      %end
      %if( ~strcmpi(atts.orf_classification, 'Verified') && ~strcmpi(atts.orf_classification, 'Uncharacterized') )
      %  if( ~strcmpi(atts.orf_classification, 'Dubious' ) )
      %    disp(['WARNING - unexpected orf_classification found - ' atts.orf_classification]);
      %  end
      %  disp('Skipping ORF characterized as ''Dubious''');
      %  return;
      %end

      if(~firstEntry)
        fprintf(results, '  ,');
      else
        fprintf(results, '  ');
      end
      attsJson = strjoin(cellfun(@(k,v) sprintf('"%s":%s',k,formatCell(v)), fieldnames(atts), struct2cell(atts), 'UniformOutput', false), ', ');

      if(isfield(atts, 'transl_table'))
        currGeneticCode = geneticcode(str2num(atts.transl_table));
        disp(['Overriding genetic code: ' atts.transl_table]);
      elseif( strcmpi(CDS.Reference, 'chrmt'))
        currGeneticCode = geneticcode(3); % Yeast mitochondrial
        disp(['Overriding genetic code: ' 3]);
      else
        currGeneticCode = geneticcode(geneticCode);
      end

      translation = nt2aa(sequence, 'GeneticCode', currGeneticCode, 'ACGTOnly', false);
      numStopCodons = sum(translation=='*');
      %disp(translation);


      % Do sanity tests
      assert(length(sequence)>3);
      %assert(mod(length(sequence),3)==0);  % TODO - This doesn't work for partial CDSs - why?
      if(~isPartialCDS)
        assert(mod(length(sequence),3)==0);
        if(~isfield(atts, 'transl_except'))
          assert(any(strcmpi(currGeneticCode.Starts, sequence(1:3)))) % Is the first a codon a start-codon?
          if(numStopCodons==0)
            disp(sprintf('Warning: no stop-codon found in entry %s', atts.ID));
            disp(sequence);
            disp(translation);
            return;  % Skip this sequence
          end
          if(numStopCodons>1)
            disp(sprintf('Error: More than one stop-codon found in entry %s', atts.ID));
            disp(sequence);
            disp(translation);
          end
          assert(numStopCodons==1); % A complete CDS must contain a stop codon (and it must be the last codon - see below)
        end
      end
      if(~isfield(atts, 'transl_except'))
        assert(numStopCodons<=1); % There is at most one stop codon

        if( numStopCodons > 0 )
          assert(translation(length(translation))=='*'); % if a stop codon exists, it must be the last codon
        end
      end



      % Write the output in json format
      fprintf(results, '{\n  "id":"%s",\n  "CDS":"%s",\n  "numExons":%d,\n  "length_nt":%d,\n  "attributes":{%s},\n  "reference":"%s",\n  "strand":"%s"\n  }\n', atts.ID, sequence, numExons, length(sequence), attsJson, CDS.Reference, CDS.Strand);
      % Also write CDS sequences in fasta format
      shortProtId = atts.ID;
      % strip the `_CDS` suffix (if present)
      suffixPos = strfind(shortProtId, '_CDS');
      if( suffixPos > 0 )
        shortProtId = shortProtId(1:suffixPos-1);
      end

      if( ~isPartialCDS)
        fastawrite(resultsFasta, shortProtId, sequence);
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
    if( ~isfield(CDS, 'dummy'))
      %if( CDS.Strand=='+')
        sequence = [sequence newseq];
      %else
      %  sequence = [newseq sequence];
      %end
    end
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
