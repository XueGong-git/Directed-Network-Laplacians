%load Dunnhumby 'The Complete Journey' data
%(https://www.dunnhumby.com/source-files/) and build the conditional
%probability matrix from the most purchased commodities
%The values represent the probability to purchase column products when row
%products are bought. If we use it as the adjecency matrix for a
%weighted directed network, the edges weights represents how likely the
%target nodes will be purchased when the source nodes is purchased.


clear all

%convert the original .csv file to .mat format. Sicne the original csv file is too big, we only uploaded the .mat format.
%filename = 'transaction_data.csv'; trans_data = readtable(filename); % Dunnhumby 'The Complete Journey' data, downloaded from https://www.dunnhumby.com/source-files/
%save('processed_data/trans_data', 'trans_data')

load('processed_data/trans_data.mat')
zeroEntry = find((trans_data.QUANTITY == 0) | (trans_data.SALES_VALUE == 0)); 
trans_data(zeroEntry,:) = [];%remove rows with zero quantity and sales value original data 

%load product information
product = readtable('datasets/product.csv');%also from Dunnhumby 'The Complete Journey' data
all_trans = join(trans_data, product);
save('processed_data/all_trans', 'all_trans')

%aggregate quantity by commodity
all_trans.C = findgroups(all_trans.COMMODITY_DESC);
aggC = splitapply(@sum,all_trans.QUANTITY,all_trans.C);

%filter basket with the high frequency commodities
[out,idx] = sort(aggC, 'descend');
nCmdty = size(aggC);
topCmdty = idx(1:ceil(10*nCmdty/100));
top_cmdty_trans = all_trans(ismember(all_trans.C, topCmdty), :);

trans = top_cmdty_trans; %use top commodities to build conditional probability matrix
minOccur = 10; %remove records of commodities that appear in less than 'minOccur' of baskets

%convert fields to group number
trans.basketID = findgroups(trans.BASKET_ID);
trans.productID = findgroups(trans.PRODUCT_ID);
trans.cmdtyID = findgroups(trans.COMMODITY_DESC);
commodity_list = unique(table(trans.cmdtyID, trans.COMMODITY_DESC));
commodity_name = unique(trans.COMMODITY_DESC);
[n_basket, ~] = size(unique(trans.basketID));
[n_product, ~] = size(unique(trans.productID));
[n_commodity, ~] = size(unique(trans.COMMODITY_DESC));

%build quantity matrix
quantity = zeros(n_basket,n_commodity);
unpivotedTdata = trans(:,{'basketID','COMMODITY_DESC','QUANTITY'});%get the unit of each commodity for each basket
pivotedTdata = unstack(unpivotedTdata, 'QUANTITY', 'COMMODITY_DESC');
%update units of products
for i = 1 : n_commodity
quantity(:, i) = splitapply(@nansum,pivotedTdata(:,i+1),pivotedTdata.basketID);
end
%convert quantity matrix to binary-e.g. a basket has A(1) or not(0)
M = (quantity>0);
colsum = sum(M, 1); %calculate the sum of each commodity across all baskets
sumVec = colsum(colsum >= minOccur); %remove records of commodities that appear in  less than 10 baskets
M = M(:,(colsum>=minOccur));
commodity_name = commodity_name((colsum>=minOccur));
n_commodity = length(sumVec);
%matrix of #(A and B) - number of occurance when A and B both appear
jointM = M.'*M;
sumM = repmat(sumVec,[length(sumVec),1]);

%conditional probability matrix P(A|B)
condP = jointM./sumM; %no of baskests where A and B appear together/ no of baskets where B appear
condP = condP'; %conP(i, j) = P(j|i)
save('datasets/condP.mat','condP', 'commodity_name');

