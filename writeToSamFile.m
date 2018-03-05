% write to fasta file
% input:
% output:
% example:
function writeToSamFile(samFormat, fileName, genomeName, genomeLength)
outfilename = fileName;
fid=fopen(outfilename,'w');
firstt = ['@SQ\tSN:', genomeName, '\t', 'LN:', num2str(genomeLength), '\n'];
fprintf(fid,firstt);
for i=1:length(samFormat)
      fprintf(fid,[samFormat{i},'\n']);
% fprintf(fid,'%s\n',Clones(i,:));
end
fclose(fid);
end