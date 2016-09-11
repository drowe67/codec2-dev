function [Hrows Hcols] = Mat2Hrows(HRA); 

H = full(HRA); 
[Nr Nc] = size(H); 

H= H(:,1:Nc-Nr); 
[Nr Nc] = size(H); 



Max_colwt = max(sum(H));
Max_rowwt = max(sum(H')); 
Hcols = zeros(Nc, Max_colwt); 
Hrows = zeros(Nr, Max_rowwt); 

for i = 1:Nr
    nz = find(H(i,:)); 
    Hrows(i,1:length(nz)) = nz; 
end

H = H'; 
for i = 1:Nc
    nz = find(H(i,:)); 
    Hcols(i,1:length(nz)) = nz; 
end
    
