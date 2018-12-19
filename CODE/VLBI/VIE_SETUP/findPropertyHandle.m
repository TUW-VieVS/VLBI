function [ hprop ] = findPropertyHandle( hobj , propname )
%FINDPROPERTYHANDLE retrieve the handle of the meta.property object

    mco = metaclass(hobj) ;
    plist = mco.PropertyList;
    for k=1:numel(plist)
        if strcmpi(plist(k).Name,propname)
%             fprintf('[%s] property was found at index #%d\n',propname,k)
            hprop = plist(k) ;
            return
        end
    end
    % no preperty was found if we are here
    hprop = [] ;
%     fprintf('[%s] property was not found.\n',propname)

end
