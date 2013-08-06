function get_table=disp_table(table_id,step_id)
% 建立資料流表格，停止於 inp, cst, var step_id

global ABST

log = ABST.log;

step_id=[step_id,0];


save_id = 0;

get_table = zeros(0,3);

get_table = data_cal(table_id,save_id,log,get_table,step_id);


function get_table = data_cal(table_id,save_id,log,get_table,step_id)

if ~isempty(find(get_table(:,1)==table_id, 1))
    return
end

id1=log(table_id,6); % 下層主
id2=log(table_id,7); % 下層輔

if log(table_id,1)==6 && log(table_id,4)==28
    % sgn, 且為參照時
    id2 = save_id;
end

ny = size(get_table,1);
get_table(ny+1,:) = [table_id id1 id2];
% get_table

% try

cond_id1=stop_cond(id1,log,step_id);
cond_id2=stop_cond(id2,log,step_id);
cond_id3=stop_cond(save_id,log,step_id);
cond = mat2str([cond_id1 cond_id2 cond_id3]);
switch cond
    case '[1 1 1]'
        return
    case '[1 1 0]'
        id1 = save_id;
        id2 = 0;
        save_id = 0;
    case '[1 0 1]'
        id1 = id2;
        id2 = 0;
        save_id = 0;
    case '[1 0 0]'
        id1 = id2;
        id2 = save_id;
        save_id = 0;
    case '[0 1 1]'
        id2 = 0;
        save_id = 0;
    case '[0 1 0]'
        id2 = save_id;
        save_id = 0;
    case '[0 0 1]'
        save_id = 0;
%     case '[0 0 0]'
end


get_table = data_cal(id1,id2,log,get_table,step_id);

if ~stop_cond(id2,log,step_id)
    get_table = data_cal(id2,save_id,log,get_table,step_id);
end

if ~stop_cond(save_id,log,step_id)
    get_table = data_cal(save_id,0,log,get_table,step_id);
end


function chk_cond=stop_cond(id,log,step_id)
% chk_cond = 1 無法採用

chk_cond = 0;
if ~isempty(find(step_id==id, 1))
    chk_cond = 1;
    return
end
switch log(id,1)
    case {1,2,5} % cst, var, inp
        chk_cond = 1;
end

