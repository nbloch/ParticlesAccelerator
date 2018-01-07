function pString = buildParamString(pNames,pSet, struct_field)

%the params string will be used as the field name in the results struct, so
%that it cannot contain any character.
% Spaces are replaced by 'X'. '.' and '-' by '_'
    if struct_field
        pSetString = strrep(string(pSet), "-", "_");
        pSetString = strrep(string(pSetString), ".", "_");
        valSeparator = 'X';
        paramSeparator = 'X';
    else
        pSetString = string(pSet);
        valSeparator = ' = ';
        paramSeparator = ', ';
    end
    pString = "";
    for i=1:numel(pSetString)
        pString = pString + pNames{i} + valSeparator + pSetString(i);
        % We don't add a separator at the end of the string
        if i < numel(pSetString)
            pString = pString + paramSeparator;
        end
    end
%     pString = strjoin(string(pNames), separator) + separator...
%         + strjoin(pSetString, separator);
end

