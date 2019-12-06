function BooleanizeFromFile( InputData, OutputFile )
%BOOLEANIZEFROMFILE Summary of this function goes here
%   Detailed explanation goes here

load('/work/users/mali/Tasks/core_networks/Evan_networks/Matlab/Thresholds_CoTFcRF.mat')
fid = fopen('/work/users/mali/Tasks/core_networks/Evan_networks/Matlab/background_genes.txt');
backgroundgenes = textscan(fid,'%s%s','CollectOutput',true);
fclose(fid);
disp(InputData);
disp(OutputFile);
fid = fopen(InputData);
inpData = textscan(fid,'%s%s%f','CollectOutput',true,'HeaderLines',1);
fclose(fid);

inputmat = -1*ones(2137,1);
for i=1:2137
    idx = find(strcmpi(backgroundgenes{1}(i),inpData{1}(:,1)));
    if(~isempty(idx))
        inputmat(i) = inpData{2}(idx);
    end
end
existing = find(inputmat > -1);
dists_inp = dists(existing);
backgroundgenes_inp = backgroundgenes{1}(existing,2);
exp_inp = inputmat(existing);
[ Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj,Pvals_inter_low_adj,Pvals_inter_up_adj ] = CalculatePValues( dists_inp, exp_inp);

s = size(exp_inp);
s = s(1);
Bool_exp = zeros(s,1);
Bool_exp(find(Pvals_low_adj > 0.1)) = 1;

output = cell(s,2);
output(:,1) = backgroundgenes_inp(:);
for i=1:s
    output{i,2} = Bool_exp(i);
end
outputTable = cell2table(output);
writetable(outputTable,OutputFile)

end

