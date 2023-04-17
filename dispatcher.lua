local num_shards = REPLACE_THIS_WITH_NUM_SHARDS
local shards_boundaries = {REPLACE_THIS_WITH_BOUNDARIES}
local shards_modem_side = "back"

local data_modem = peripheral.wrap(shards_modem_side)

term.clear()
term.setCursorPos(1, 1)

local function wait_for_shards()
    print("Dispatcher startup procedure")

    for i=1, num_shards do
        print("Waiting for", i, "shard.")
        data_modem.open(i)

        while true do
            local _, side, channel, _, message, _ = os.pullEvent("modem_message")

            if side == shards_modem_side and channel == i then
                local msg = textutils.unserialize(message)
                if msg["is_startup_procedure"]
                        and msg["is_shard"]
                        and msg["num_shard"] == i  then
                    local response = {is_startup_procedure=true, is_shard=false, num_shard=i}
                    print("Got response from", i, "shard")

                    data_modem.transmit(i, i, textutils.serialize(response))
                    break;
                end
            end
        end

        data_modem.close(i)
    end
end

local function get_shard(y)
    for k, v in pairs(shards_boundaries) do
        -- 1 - min_y incl, 2 - max_y not incl
        if y >= tonumber(v[1]) and y < tonumber(v[2]) then
            return k
        end
    end
    return -1
end

local function wait_for_reply(shard_num, x, y)
    while true do
        local _, side, channel, _, message, _ = os.pullEvent("modem_message")
        local msg = textutils.unserialize(message)
        if side == shards_modem_side and channel == shard_num
                and msg["is_response"]
                and msg["num_shard"] == shard_num
                and msg["y_pos"] == y
                and msg["x_pos"] == x
        then
            msg["is_response"] = nil
            msg["num_shard"] = nil
            return msg
        end
        os.sleep(0)
    end
end

local function request_data(request, channel)
    while true do
        data_modem.transmit(channel, channel, textutils.serialize(request))
        os.sleep(0)
    end
end

local function get_data(x, y)
    local shard_num = get_shard(y)

    if shard_num <= 0 then
        return {
            x_pos=x,
            y_pos=y,
            delta_t=-1,
            pitch=-1,
            airtime=-1
        }
    end

    local request = {
        num_shard=shard_num,
        y_pos=y,
        x_pos=x,
        is_request=true
    }

    data_modem.open(shard_num)

    local response = {}

    parallel.waitForAny(
            (function () request_data(request, shard_num) end),
            (function () response = wait_for_reply(shard_num, x, y) end)
    )

    data_modem.close(shard_num)

    return response
end

local function wait_for_request()
    while true do
        print("write y")
        local y = tonumber(io.read())
        print("write x")
        local x = tonumber(io.read())

        print(textutils.serialize(get_data(x, y)))
    end
end

wait_for_shards()
wait_for_request()
