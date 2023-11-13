local data_modem = peripheral.find("modem")
data_modem.open(REPLACE_THIS_WITH_MODEM_CHANNEL)

local function check_if_header(name)
    local start, stop = string.find(name, "REPLACE_THIS_WITH_HEADER_NAME")
    if start ~= 1 then return false end
    return true
end

local headers = {}
for k, filename in pairs(fs.list("/")) do
    if check_if_header(filename) then
        table.insert(headers, filename)
    end
end

if #headers == 0 then error("No data file detected")
else
    term.clear()
    term.setCursorPos(1, 1)
    print("Data loaded:")
    for i=1, #headers do
        term.setCursorPos(1, i+1)
        print(i..") "..headers[i])
    end
end

while true do
while true do
    local _, _, _, _, req, _ = os.pullEvent("modem_message")
    if req == nil then break end
    if req.action == nil then break end
    if req.action == "give_headers" then
        for k, header in pairs(headers) do
            data_modem.transmit(REPLACE_THIS_WITH_MODEM_CHANNEL, REPLACE_THIS_WITH_MODEM_CHANNEL, {["header"]=header})
        end
    elseif req.action == "give_data" then
        for k, header in pairs(headers) do
            if header == req.header then
                local h = fs.open(header, "r")
                data_modem.transmit(REPLACE_THIS_WITH_MODEM_CHANNEL, REPLACE_THIS_WITH_MODEM_CHANNEL, {["header"]=header, data=h.readAll()})
                h.close()
            end
        end
    end
end
end