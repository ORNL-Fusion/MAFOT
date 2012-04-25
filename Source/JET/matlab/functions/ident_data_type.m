function [Laminar, Target] = ident_data_type(FileToOpen,Laminar, Target)
%[Laminar, Target] = ident_data_type(FileToOpen,Laminar, Target)
%Funktion read filename and try to identify the data type
%Laminar, Target = -1 activates automatic mode, respectively

    if (Laminar == -1)
        if (findstr(FileToOpen,'lam_')~=0), Laminar = 1;
        elseif (findstr(FileToOpen,'foot_')~=0), Laminar = 0;
        else
            error('ident_data_type:type_not_found','ERROR in automatic mode: File cannot be identified (Laminar or not?!).\nUse manual mode!!!');
        end
    end
    if ((Laminar == 0) && (Target == -1))
        if (findstr(FileToOpen,'_inn_')~=0), Target = 0;
        elseif (findstr(FileToOpen,'_out_')~=0), Target = 1;
        else
            error('ident_data_type:type_not_found','ERROR in automatic mode: File cannot be identified (inner or outer target?!).\nUse manual mode!!!');
        end
    end 
end