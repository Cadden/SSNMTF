function write_infile( filename,source_matrix )
%  write matrix to file
fid=fopen(filename,'w');

[x,y]=size(source_matrix);

for i=1:x
    for j=1:y
        if source_matrix(i,j) == 0
            break;
        else
            fprintf(fid,'%d\t',source_matrix(i,j));
        end
    end
    fprintf(fid,'%d\n', '');%每一行回车\n
end

fclose(fid);

end
