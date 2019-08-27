function S = read_benson(name)
    base_path = '/home/govindjsk/repos/data/benson/hypergraphs/';
    path = strcat(base_path, name, '/', name, '.mat');
    S = load(path, 'S');
end