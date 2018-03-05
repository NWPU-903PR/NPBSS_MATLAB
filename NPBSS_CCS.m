function NPBSS_CCS(referenceFile,ccsfastqFile)
if nargin < 2
    fprintf('\nExample: NSSPB(''genome.fa'', ''ccs.fastq'')\n');
    fprintf('There need at least 2 input arguments.\nOne is reference file, the other is CCS fastq file.\n');
    %fprintf('-ins(optional)          Insertion error rate between (0,1),(default: 0.00939).\n');
    %fprintf('-del(optional)          Deletion error rate between (0,1),(default: 0.00957).\n');
    %fprintf('-sub(optional)          Substitution error rate between (0,1),(default: 0.0374).\n');
    return;
end
insertMean = 0.00939;
deleteMean = 0.00957;
subsituteMean = 0.0374;
fileName = referenceFile;

[head, genome] = fastaread(referenceFile);

[heads, seqis, qvss] = fastqread(ccsfastqFile);
seq_num = length(seqis);
list = zeros(1, seq_num);
for i = 1 : seq_num
    list(i) = length(seqis{i});
end
%load('CCS_length.mat')
%load('CCS_qvis')
% sample_rate = genomeLength * depth/sum(seq_length);
% sample_num = floor(sample_rate * length(seq_length));
position = 1 : seq_num;
[rightSeq, simulatedSeq, qv] = simulatePacBio_CCS(list, genome, position, qvss, insertMean, deleteMean, subsituteMean);
writeToFastaFile(simulatedSeq,[fileName,'.npbss_simulated_CCS.fa']);
writeToFastaFile(rightSeq,[fileName,'.reads_correct_CCS.fa']);
writeToFastqFile(simulatedSeq, qv, [fileName, '.npbss_simulated_CCS.fq']);
end