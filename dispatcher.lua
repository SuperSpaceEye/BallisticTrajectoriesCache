local needed_headers = REPLACE_THIS_WITH_NEEDED_HEADERS
local data_modem_side = "back"
local data_modem = peripheral.wrap(data_modem_side)
data_modem.open(REPLACE_THIS_WITH_MODEM_CHANNEL)

local function nKeys(t) local m=0 for _,_ in pairs(t) do m=m+1 end return m end

local function header_to_offset(name)
    local start, stop = string.find(name, "REPLACE_THIS_WITH_HEADER_NAME")
    if start ~= 1 then return false, nil end
    local header = string.sub(name, stop+1)
    return true, header
end

local function await_all_headers()
    local headers_yet_to_verify = {}
    local got_data = {}
    for k, v in pairs(needed_headers) do headers_yet_to_verify["REPLACE_THIS_WITH_HEADER_NAME"..v]=true end

    parallel.waitForAll(
    function()
       while nKeys(headers_yet_to_verify) ~= 0 do
           data_modem.transmit(REPLACE_THIS_WITH_MODEM_CHANNEL, REPLACE_THIS_WITH_MODEM_CHANNEL, {action="give_headers"})
           sleep()
       end
    end,
    function()
       while nKeys(headers_yet_to_verify) ~= 0 do
       while nKeys(headers_yet_to_verify) ~= 0 do
           local _, _, _, _, req, _ = os.pullEvent("modem_message")
           if req == nil then break end
           if req.header == nil then break end
           headers_yet_to_verify[req.header] = nil
       end
       end
    end
    )
    for k, header in pairs(needed_headers) do
        header = "REPLACE_THIS_WITH_HEADER_NAME"..header
        while true do
            data_modem.transmit(REPLACE_THIS_WITH_MODEM_CHANNEL, REPLACE_THIS_WITH_MODEM_CHANNEL, {action="give_data", ["header"]=header})
            local _, _, _, _, req, _ = os.pullEvent("modem_message")
            if req ~= nil and req.header == header and req.data ~= nil then
                local converted, offset = header_to_offset(header)
                if converted then got_data[offset] = req.data end
                break
            end
        end
    end
    return got_data
end

local function combine_data(data)
    local combined_data = {}
    for k, part in pairs(data) do
        local start = tonumber(k)
        part = textutils.unserialise(part)
        print(start, nKeys(part))
        --if REPLACE_THIS_WITH_IF_RECOMBINE then
        --    for i, line in pairs(part) do
        --        local x = 2
        --        local newline = {line[1]}
        --
        --        while x <= #line do
        --            table.insert(newline, {line[x], line[x+1], line[x+2]})
        --            x = x + 3
        --        end
        --        part[i] = newline
        --    end
        --end

        if REPLACE_THIS_WITH_IF_RECOMBINE then
            for i, line in ipairs(part) do
                local x = 2
                local nx = line[1]
                if i - 256 + start == 0 then
                    print(i-256+start, i+start, i, start, nx)
                end
                local newline = {}

                while x <= #line do
                    newline[nx] = {line[x], line[x+1], line[x+2]}
                    x = x + 3
                    nx = nx + 1
                end
                part[i] = newline
            end
        end

        for i=start, nKeys(part)+start-1 do
            if part[i-start+1] == nil then print("A FUCKING ERROR"); error() end
            combined_data[i] = part[i-start+1]
        end
    end
    return combined_data
end

local data = combine_data(await_all_headers())

local get_data = REPLACE_THIS_WITH_IF_RECOMBINE and (
        function(y, x)
            local line = data[y+REPLACE_THIS_WITH_Y_DISPLACEMENT]
            --return line[x-line[1]] --get x displacement from first line value, do math
            return line[x] --get x displacement from first line value, do math
        end) or (
        function (y, x)
            local line = data[y+REPLACE_THIS_WITH_Y_DISPLACEMENT]
            x = (x-line[1])*3 + 2 --get x displacement from first line value, do math
            return {line[x], line[x+1], line[x+2]}
        end)

print(REPLACE_THIS_WITH_Y_DISPLACEMENT)
print(textutils.serialize(get_data(0, 38)))

return {
    data = data,
    get_data = get_data,
    y_displacement = REPLACE_THIS_WITH_Y_DISPLACEMENT,
    recombine_at_runtime = REPLACE_THIS_WITH_IF_RECOMBINE
}