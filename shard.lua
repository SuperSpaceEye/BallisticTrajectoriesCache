local data = {REPLACE_THIS_WITH_DATA}
local shard_num = REPLACE_THIS_WITH_SHARD_NUM
local displacement = REPLACE_THIS_WITH_SHARD_DISPLACEMENT
local shard_modem_side = "back"

local data_modem = peripheral.wrap(shard_modem_side)
data_modem.open(shard_num)

term.clear()
term.setCursorPos(1, 1)

local function push_init()
    while true do
        local message = {is_startup_procedure=true, is_shard=true, num_shard=shard_num}
        data_modem.transmit(shard_num, shard_num, textutils.serialize(message))
        os.sleep(0)
    end
end

local function pull_init()
    while true do
        local _, side, channel, _, message, _ = os.pullEvent("modem_message")
        if side == shard_modem_side and channel == shard_num then
            local message = textutils.unserialize(message)
            if message["is_startup_procedure"]
                    and not message["is_shard"]
                    and message["num_shard"] == shard_num then
                print("Got response from dispatcher.")
                return;
            end
        end
        os.sleep(0)
    end
end

local function wait_for_dispatcher()
    print("Shard", shard_num ,"startup procedure")
    parallel.waitForAny(pull_init, push_init)
end

-- Request message  {num_shard, y_pos, x_pos, is_request}
-- Response message {num_shard, y_pos, x_pos, delta_t, pitch, airtime, is_response}

local function wait_for_request()
    while true do
        local _, side, channel, _, message, _ = os.pullEvent("modem_message")
        print(message)
        if side == shard_modem_side and channel == shard_num then
            local msg = textutils.unserialize(message)
            if msg["is_request"] and msg["num_shard"] == shard_num then
                local response = {num_shard=shard_num,
                                  y_pos=msg["y_pos"],
                                  x_pos=msg["x_pos"],
                                  delta_t=-1, pitch=-1, airtime=-1,
                                  is_response=true
                }

                local line = data[tonumber(msg["y_pos"])+displacement]

                if line ~= nil then
                    local item = line[tonumber(msg["x_pos"])]
                    if item ~= nil then
                        response = {num_shard=shard_num,
                                    y_pos=msg["y_pos"],
                                    x_pos=msg["x_pos"],
                                    delta_t=item[1], pitch=item[2], airtime=item[3],
                                    is_response=true
                        }
                    end
                end

                data_modem.transmit(shard_num, shard_num, textutils.serialize(response))
            end
        end
    end
end

wait_for_dispatcher()
wait_for_request()
