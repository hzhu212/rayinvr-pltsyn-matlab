function [ws] = run2struct(script)
% 运行一个脚本，并将其局部变量保存到一个 struct 中

    run(script);
    ws = ws2struct();
end


function [ws] = ws2struct()
% 将 workspace 转化为 struct

    ws = struct();
    vars = evalin('caller', 'who');
    for ii = 1:numel(vars)
        ws.(vars{ii}) = evalin('caller', vars{ii});
    end
end
