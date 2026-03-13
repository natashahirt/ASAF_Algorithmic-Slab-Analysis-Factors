# Global reference to keep the server alive
const server_ref = Ref{Any}(nothing)

function start_websocket()
    if !isnothing(server_ref[])
        println("Server already running")
        return
    end
    
    server_ref[] = WebSockets.listen!("127.0.0.1", 2000) do ws
        for msg in ws
            if msg == "init"
                continue
            end
            line_coordinates = get_tributary_lines(msg)
            WebSockets.send(ws, line_coordinates)
            println("Sent")
        end
    end
    
    println("WebSocket server started on ws://127.0.0.1:2000")
end

function stop_websocket()
    if !isnothing(server_ref[])
        close(server_ref[])
        server_ref[] = nothing
        println("WebSocket server stopped")
    end
end

function restart_websocket()
    stop_websocket()
    start_websocket()
end