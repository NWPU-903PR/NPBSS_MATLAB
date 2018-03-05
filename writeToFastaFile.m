% write to fasta file
% input:
% output:
% example:
function writeToFastaFile(simulatedSeqs, fileName)
outfilename = fileName;
fid=fopen(outfilename,'w');

for i=1:length(simulatedSeqs)
      fprintf(fid,'>%s\n%s\n',['Simulated_' num2str(i)],simulatedSeqs{i});
% fprintf(fid,'%s\n',Clones(i,:));
end
fclose(fid);
end