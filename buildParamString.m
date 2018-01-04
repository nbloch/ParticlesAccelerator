function pString = buildParamString(pNames,pSet, struct_field)

%the params string will be used as the field name in the results struct, so
%that it cannot contain any character.
% Spaces are replaced by 'X'. '.' and '-' by '_'
    if struct_field
        pSetString = strrep(string(pSet), "-", "_");
        pSetString = strrep(string(pSetString), ".", "_");
        separator = 'X';
    else
        pSetString = string(pSet);
        separator = ' ';
    end
    
    pString = strjoin(string(pNames), separator) + separator...
        + strjoin(pSetString, separator);
end

