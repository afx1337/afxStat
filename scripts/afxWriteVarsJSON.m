function afxWriteVarsJSON(fname, varargin)
    data = struct();

    for k = 1:length(varargin)
        name = inputname(k+1);

        if isempty(name)
            name = ['var' num2str(k)];
        end

        data.(name) = varargin{k};
    end

    jsonStr = prettyjson(jsonencode(data));
    
    fid = fopen(fname,'w');
    if fid == -1, error('File not writeable.'); end
    fprintf(fid,'%s',jsonStr);
    fclose(fid);
end

function pretty = prettyjson(str)

    indentSize = 2;
    indent = 0;
    pretty = '';
    inString = false;

    for i = 1:length(str)
        ch = str(i);

        if ch == '"' && (i == 1 || str(i-1) ~= '\')
            inString = ~inString;
        end

        if ~inString
            switch ch
                case {'{','['}
                    indent = indent + 1;
                    pretty = [pretty ch newline repmat(' ',1,indent*indentSize)];
                case {'}',']'}
                    indent = indent - 1;
                    pretty = [pretty newline repmat(' ',1,indent*indentSize) ch];
                case ','
                    pretty = [pretty ch newline repmat(' ',1,indent*indentSize)];
                case ':'
                    pretty = [pretty ch ' '];
                otherwise
                    pretty = [pretty ch];
            end
        else
            pretty = [pretty ch];
        end
    end
end