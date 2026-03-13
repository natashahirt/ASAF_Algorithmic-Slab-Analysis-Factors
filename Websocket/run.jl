# Globals to manage the WebSocket server
const websocket_server_task = Ref{Union{Nothing, Task}}(nothing)
const server_running = Ref(false)  # A flag to indicate server state
const websocket_server = Ref{Union{Nothing, WebSockets.WebSocket}}(nothing)
const active_connections = Ref{Set{WebSockets.WebSocket}}(Set{WebSockets.WebSocket}())
const message_queue = Ref{Vector{Dict{String,Any}}}(Vector{Dict{String,Any}}())

"""
    start_websocket_server(host::String, port::Int)

Starts a WebSocket server on the specified host and port.

Args:
    host: The host address to bind the server to.
    port: The port number to listen on.
"""
function start_websocket_server(; host::String = "127.0.0.1", port::Int = 2000)
    if server_running[]
        # Check if the server task is actually running
        if !isnothing(websocket_server_task[]) && !istaskdone(websocket_server_task[])
            println("Server is already running.")
            return
        else
            # Reset the state if the task is done but flag wasn't updated
            server_running[] = false
            websocket_server[] = nothing
            websocket_server_task[] = nothing
            println("Detected stale server state, resetting...")
        end
    end

    websocket_server_task[] = @async begin
        server_running[] = true
        try
            println("Starting WebSocket server on ws://$host:$port...")
            websocket_server[] = WebSockets.listen!(host, port) do ws
                try
                    println("New WebSocket connection established.")
                    push!(active_connections[], ws)
                    
                    while isopen(ws) && server_running[]
                        msg = WebSockets.receive(ws)
                        if msg == "init"
                            println("Received init message")
                            continue
                        end
                        msg_dict = JSON.parse(msg)
                        push!(message_queue[], msg_dict)
                        println("Received message: ", msg_dict)
                        WebSockets.send(ws, "Received")
                    end
                catch e
                    if !server_running[]
                        println("Connection closed during shutdown")
                    else
                        println("Error handling WebSocket connection: ", e)
                    end
                finally
                    delete!(active_connections[], ws)
                end
            end
            wait(websocket_server[])
        catch e
            if isa(e, Base.IOError) && occursin("address already in use", sprint(showerror, e))
                println("Port $port is already in use. Try a different port.")
            else
                println("Failed to start WebSocket server: ", e)
                println("Error details:")
                showerror(stdout, e)
                println()
            end
        finally
            server_running[] = false
            println("Server has stopped.")
        end
    end

    println("WebSocket server started on ws://$host:$port")
end

"""
    stop_websocket_server()

Stops the running WebSocket server.
"""
function stop_websocket_server()
    if !isnothing(websocket_server[]) || server_running[]
        server_running[] = false
        
        # Close all active connections first
        for ws in active_connections[]
            try
                if isopen(ws)
                    close(ws)
                end
            catch e
                @warn "Error closing websocket connection: $e"
            end
        end
        empty!(active_connections[])
        
        # Force close the HTTP server and its TCP listener
        if !isnothing(websocket_server[])
            try
                # Get the underlying TCP listener
                listener = WebSockets.getsocket(websocket_server[])
                # Close the HTTP server
                close(websocket_server[])
                # Explicitly close the TCP listener
                close(listener)
                # Force GC to clean up the closed connections
                GC.gc()
            catch e
                @warn "Error closing websocket server: $e"
            end
        end

        # Attempt to interrupt the server task if it exists
        if !isnothing(websocket_server_task[]) && !istaskdone(websocket_server_task[])
            try
                schedule(websocket_server_task[], InterruptException(), error=true)
                # Give it a short time to clean up
                sleep(0.5)
            catch e
                @warn "Error interrupting server task: $e"
            end
        end

        # Reset all state variables
        websocket_server[] = nothing
        websocket_server_task[] = nothing
        
        println("WebSocket server stopped.")
    else
        println("No WebSocket server running.")
    end
end
