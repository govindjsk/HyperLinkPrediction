function [S_train, S_test, non_Ss] = prepare_benson_data(data_name)
    path = strcat('/home/govindjsk/repos/datagen/', data_name, '.mat');
    S = load(path);
    Strain = sparse(double(S.Strain));
    Stest = sparse(double(S.Stest));
    Negative = S.Negative;
    clear S;
    for i = 1:length(Negative)
        Negative{i} = sparse(double(Negative{i}));
        Negative{i}(Negative{i} > 0) = 1;
    end
    S_train = Strain;
    S_test = Stest;
    non_Ss= Negative;
end

% function B = convert_to_sparse(A)
%     [m, n] = size(A);
%     nz = A~=0;
%     [rows, columns] = find(nz);
%     values = A(nz);
%     B = sparse(rows, columns, values, m, n);
% end