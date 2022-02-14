function letter=alphabet(index,capital)

if ~exist('capital','var')
    capital=false(size(index));
end
list='abcdefghijklmnopqrstuvwxyz';

letter=list(index);

letter(capital)=upper(letter(capital));
    