% write to fastq file
% input:
% output:
% example:
function writeToFastqFile(simulatedSeq, qv, fileName)
outfilename = fileName;
fid=fopen(outfilename,'w');

for i=1:length(simulatedSeq)
    q = qv{i};
    phred = q + 33;
    strr = char(phred);
    fprintf(fid,'@%s\n%s\n%s\n%s\n',['Simulated_' num2str(i)],simulatedSeq{i},'+',strr);
% fprintf(fid,'%s\n',Clones(i,:));
end
fclose(fid);
end