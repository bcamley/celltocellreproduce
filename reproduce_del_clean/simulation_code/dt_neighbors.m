function neighs = dt_neighbors(DT,N)
neighs = cell(1,N);
for i = 1:N
    neighs{i} = [];
end

ed = DT.edges;
for k = 1:size(ed,1)
    neighs{ed(k,1)} = [neighs{ed(k,1)} ed(k,2)];
    neighs{ed(k,2)} = [neighs{ed(k,2)} ed(k,1)];
end