%  Main program of hyperlink prediction experiments on metabolic networks
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
clearvars;

% 1: email-Enron, 2: NDC-substances, 3: contact-primary-school, 4:
% tags-math-sx, 5: DBLP
data_names = {'email-Enron', 'NDC-substances', 'contact-primary-school', ...
    'tags-math-sx', 'DBLP'};
data_ids = 1:length(data_names);
data_id_name_map = containers.Map(data_ids, data_names);

% Method--0: SCCMM  1: CMM  2: BS  3: SHC  4: HPLSF  5: FM  6: Katz  7: CN 
%                  (8: LR  9: NN  10: DPP, not included in paper)
method_names = {'SCCMM', 'CMM', 'BS', 'SHC', 'HPLSF', 'FM', 'Katz', 'CN', ...
    'LR', 'NN', 'DPP'};
method_ids = 0:length(method_names)-1;
method_id_name_map = containers.Map(method_ids, method_names);

data_ids = [1];

for data_id = data_ids
    data_name = data_id_name_map(data_id)
    [S_train, S_test, non_Ss] = prepare_benson_data(data_name);
    num_missing = size(S_test, 2);
    num_exp = 1;
    num_neg_samplings= min([1, length(non_Ss)]);
    method_ids = [7]; % modify this list to compare more methods
    for m = 1:length(method_ids)
        method_id = method_ids(m);
        method_name = method_id_name_map(method_id)
        matched_hl_count = zeros(num_neg_samplings, num_exp);
        selected_hl_count = zeros(num_neg_samplings, num_exp);
        matched_hl_avg = zeros(num_neg_samplings, num_exp);
        hl_pred_auc = zeros(num_neg_samplings, num_exp);
        for t = 1:num_neg_samplings
            non_S = non_Ss{t};
            %poolobj = parpool(feature('numcores')); % uncomment this line and change the "for" in next line to "parfor" in order to run in parallel
            for e = 1:num_exp
                e
                train = S_train;
                num_train = size(train, 2);
                test = [S_test, non_S];
                num_test = size(test, 2);
                train_ids = 1:num_train;
                test_ids = 1:num_test;
                pos_test_ids = 1:num_missing;
                neg_test_ids = (num_missing+1): num_test;
                labels = [ones(1, size(S_test, 2)), zeros(1, size(non_S, 2))];                
                
                % METHOD
                % Lambda: the 1/0 indicator vector indicating if each test hyperlink is predicted positive
                % scores: the soft indicator vector recording the confidence score of each test hyperlink being positive
%                 method_name = method_id_name_map(method_id);
                
                [Lambda,scores] = HLpredict(train, test, num_missing,method_name,e,labels);
                predictions = test_ids(Lambda');
                true_pos = intersect(pos_test_ids, predictions);
%                 predictions = test_names(Lambda');
%                 match = strcmp(repmat(missing_names,1,size(predictions,2)),repmat(predictions,size(missing_names,1),1)); 
%                 num_of_matched_reactions = nnz(match) % see how many predictions are real missing reactions
                true_pos_count = length(true_pos);
                matched_hl_count(t,e) = true_pos_count;

%                 assert(nnz(Lambda)==tau)

                num_of_selected_reactions = nnz(Lambda);
                selected_hl_count(t,e) = num_of_selected_reactions;
                average_guess_match_num = num_missing*num_missing/size(Lambda,1);
                matched_hl_avg(t,e) = average_guess_match_num;
                
                % calculate AUC
                [~,~,~,auc] = perfcurve(labels,scores,true);
                hl_pred_auc(t,e) = auc;
            end
            if exist('poolobj')
                delete(poolobj)
            end
            
        end

        % calculate mean results
        Recall = bsxfun(@times,matched_hl_count,1./num_neg_samplings');
        Precision = matched_hl_count./selected_hl_count;
        average_recall = mean(Recall,2);
        std_recall = std(Recall,0,2);
        average_precision = mean(Precision,2);
        std_precision = std(Precision,0,2);
        average_guess_match_num = mean(matched_hl_avg,2)
        std_guess_match_num = std(matched_hl_avg,0,2);
        average_AUC = mean(hl_pred_auc,2)
        std_AUC = std(hl_pred_auc,0,2);
        average_match_num = mean(matched_hl_count,2)
        std_match_num = std(matched_hl_count,0,2);
        sound(sin(2*pi*25*(1:4000)/100));  % make a sound in the end
        %% save results
        save(sprintf('result/%s_%s.mat',data_name,method_name),'average_match_num', 'std_match_num', 'average_guess_match_num', 'std_guess_match_num', 'average_AUC', 'std_AUC','average_recall','std_recall','average_precision','std_precision','num_missing');
    end
end
