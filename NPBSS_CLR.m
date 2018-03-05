function NPBSS_CLR(referenceFile,pnames)
if nargin < 2
    fprintf('-----------------------------------------\n');
    fprintf('NPBSS needs at least 2 input arguments\n');
    fprintf('Example: NPBSS_CLR(''EcoliGenome.fa'', ''-dep 5'')\n');
    fprintf('---Usage:\n');
    fprintf('-dep          Sequencing depth.\n');
    fprintf('-n            Number of sequences.\n');
    fprintf('-min          Minimum sequence length (default:100).\n');
    fprintf('-max          Maximum sequence length (default:45000).\n');
    fprintf('-samp         Sample the read length from a FASTA (Ex.: -samp fineName.fa).\n');
    fprintf('-lg mean std  The mean and standard deviation for log-normal (Ex.: -lg 8500 6930).\n');
    fprintf('-len          Average sequence length (default: 8500).\n');
    fprintf('-sub          Substitution error rate between (0,1),(default: 0.06).\n');
    fprintf('-ins          Insertion error rate between (0,1),(default: 0.03).\n');
    fprintf('-sam          Sam format output (default: 0).\n');
    fprintf('-del          Deletion error rate between (0,1),(default: 0.06).\n');
    fprintf('-qv           QVs selection table, not recommended to change.\n');
    fprintf('-model        Error model offered by users.\n');
    return;
end

s = regexp(pnames, '\s+', 'split');

depth = 0;
n = 0;
meanis = 8500;
varis  = 6953;
ins_is = 0.02953;
del_is = 0.06403;
sub_is = 0.0230;
samFlag = 0;
minLen = 100;
maxLen = 45000;
qvNumPorpation = [0.000276726086402811,0.00832439676688714,0.0253280741491020,0.0323037535959514,0.0306580764541609,...
                  0.0359438422593862,0.0405823269033696,0.0479061456202096,0.0606710538165784,0.0689325673079174,...
                  0.0724403432117603,0.0895914094701023,0.121460336125986,0.206006999285471,0.159571766877217];
model = [0.3634,0.3421,0.2868,0.1766,0.1482,0.1127,0.0978,0.0869,0.0686,0.0584,0.0555,0.0486,0.0426,0.0307,0.0207];
for i = 1:length(s)
    aa = s{i};
    if strcmp(aa, '-dep')
        depth = str2num(s{i + 1});
        %break;
    elseif strcmp(aa, '-lg')
        cc = s{i + 1};
        ccc = regexp(cc, '\s+', 'split');
        meanis = ccc{1};
        varis  = ccc{2};
        %break
    elseif strcmp(aa, '-n')
        n = str2double(s{i + 1});
    elseif strcmp(aa, '-samp')
        sampe_file = s{i + 1};
    elseif strcmp(aa, '-qv')
        qvNumPorpation = str2num(s{i + 1});
    elseif strcmp(aa, '-model')
        model = str2num(s{i + 1});
    elseif strcmp(aa, '-len')
        meanis = str2double(s{i + 1});
    elseif strcmp(aa, '-min')
        minLen = str2double(s{i + 1});
    elseif strcmp(aa, '-max')
        maxLen = str2double(s{i + 1});
    elseif strcmp(aa, '-sub')
        sub_is = str2double(s{i + 1});
    elseif strcmp(aa, '-ins')
        ins_is = str2double(s{i + 1});
    elseif strcmp(aa, '-del')
        del_is = str2double(s{i + 1});
    elseif strcmp(aa, '-sam')
        samFlag = 1;
    end
end

if depth <= 0 && n <=0
    if meanis > 0 && varis > 0
        depth = 5;
    elseif n > 0
        [head, genome] = fastaread(referenceFile);
        depth = ceil(8500 * n/length(genome));
    elseif ischar(sampe_file)
        [head, seq] = fastaread(sampe_file);
        n = length(seq);
        lenIs = zeros(n,1);
        for i = 1:n
            lenIs(i) = length(seq{i});
        end
        meanis = mean(lenIs);
        varis  = std(lenIs);
        [head, genome] = fastaread(referenceFile);
        depth = ceil(meanis * n/length(genome));
    end
end   


fileName = referenceFile;
[head, genome] = fastaread(referenceFile);
genomeName = head;
if n > 0
    nSimulated = n;
else 
    nSimulated = floor(depth * length(genome)/meanis);
end
%genomeLength = length(genome);
%load('ecoilLength.mat')
m = meanis;
v = varis^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
list1 = ceil(lognrnd(mu,sigma,1,nSimulated));
list = list1(list1 >= minLen);
list = list(list <= maxLen);
%list = generateLengthList(ecoilLength, genomeLength, depth);
[rightSeq, simulatedSeq, qv, samFormat] = simulatePacBio4(list, genome, genomeName, ins_is, del_is, sub_is, samFlag, qvNumPorpation, model);
writeToFastaFile(simulatedSeq,[fileName,'.npbss_simulated_CLR.fa']);
writeToFastaFile(rightSeq,[fileName,'.reads_correct_CLR.fa']);
writeToFastqFile(simulatedSeq, qv, [fileName, '.npbss_simulated_CLR.fq']);
if samFlag
    writeToSamFile(samFormat, [fileName,'.npbss_simulated_CLR.sam'], genomeName, length(genome));
end
end