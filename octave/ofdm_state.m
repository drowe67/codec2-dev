% ofdm_state.m
%
% Library of state machine functions for the OFDM modem

1;

%-------------------------------------------------------------------
% sync_state_machine - calls mode-specific sync state state_machine
%-------------------------------------------------------------------

function states = sync_state_machine(states, rx_uw)
  if strcmp(states.state_machine, "voice1")
    states = sync_state_machine_voice1(states, rx_uw);
  elseif strcmp(states.state_machine, "data")
    if strcmp(states.data_mode, "streaming")
      states = sync_state_machine_data_streaming(states, rx_uw);
    else
      states = sync_state_machine_data_burst(states, rx_uw);
    end
  elseif strcmp(states.state_machine, "voice2")
    states = sync_state_machine_voice2(states, rx_uw);
  else
    assert(0);
  endif
endfunction

%--------------------------------------------------------------------
%  Due to the low pilot symbol insertion rate and acquisition issues
%  the earlier OFDM modem waveforms (700D and 2020) need a complex
%  state machine to help them avoid false sync.
%--------------------------------------------------------------------

function states = sync_state_machine_voice1(states, rx_uw)
  ofdm_load_const;
  next_state = states.sync_state;
  states.sync_start = states.sync_end = 0;

  if strcmp(states.sync_state,'search')

    if states.timing_valid
      states.frame_count = 0;
      states.sync_counter = 0;
      states.modem_frame = 0;
      states.sync_start = 1;
      next_state = 'trial';
    end
  end

  if strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial')

    states.frame_count++;

    % UW occurs at the start of a packet
    if states.modem_frame == 0
        states.uw_errors = sum(xor(tx_uw,rx_uw));

        if strcmp(states.sync_state,'trial')
          if states.uw_errors >= states.bad_uw_errors
            states.sync_counter++;
            states.frame_count = 0;
          end
          if states.sync_counter == 2
            next_state = "search";
            states.phase_est_bandwidth = "high";
          end
          if states.frame_count == 4
            next_state = "synced";
            % change to low bandwidth, but more accurate phase estimation
            states.phase_est_bandwidth = "low";
          end
          if states.uw_errors < 2
            next_state = "synced";
            % change to low bandwidth, but more accurate phase estimation
            states.phase_est_bandwidth = "low";
          else
            next_state = "search";
          end
        end

        if strcmp(states.sync_state,'synced')
          if states.uw_errors > 2
            states.sync_counter++;
          else
            states.sync_counter = 0;
          end

          if states.sync_counter == 6
            next_state = "search";
            states.phase_est_bandwidth = "high";
          end
        end
      end % if modem_frame == 0 ....

      % keep track of where we are up to in packet
      states.modem_frame++;
      if (states.modem_frame >= states.Np) states.modem_frame = 0; end
  end

  states.last_sync_state = states.sync_state;
  states.sync_state = next_state;
endfunction


%-------------------------------------------------------
% data (streaming mode) state machine
%-------------------------------------------------------

function states = sync_state_machine_data_streaming(states, rx_uw)
  ofdm_load_const;
  next_state = states.sync_state;
  states.sync_start = states.sync_end = 0;

  if strcmp(states.sync_state,'search')
    if states.timing_valid
      states.sync_start = 1; 
      states.sync_counter = 0;
      next_state = 'trial';
    end
  end

  states.uw_errors = sum(xor(tx_uw,rx_uw));

  if strcmp(states.sync_state,'trial')
    if states.uw_errors < states.bad_uw_errors;
      next_state = "synced";
      states.packet_count = 0;
      states.modem_frame = Nuwframes;
    else
      states.sync_counter++;
      if states.sync_counter > Np
        next_state = "search";
      end
    end
  end
 
  % Note packetsperburst==0 we don't ever lose sync, which is useful for 
  % stream based testing or external control of state machine
  
  if strcmp(states.sync_state,'synced')
    states.modem_frame++;
    if (states.modem_frame >= states.Np) 
      states.modem_frame = 0; 
      states.packet_count++;
      if (states.packetsperburst)
        if (states.packet_count >= states.packetsperburst)
          next_state = "search";
        end
      end
    end
  end
  
  states.last_sync_state = states.sync_state;
  states.sync_state = next_state;
endfunction

%-------------------------------------------------------
% data (burst mode) state machine
%-------------------------------------------------------

function states = sync_state_machine_data_burst(states, rx_uw)
  ofdm_load_const;
  next_state = states.sync_state;
  states.sync_start = states.sync_end = 0;

  if strcmp(states.sync_state,'search')
    if states.timing_valid
      states.sync_start = 1; 
      states.sync_counter = 0;
      next_state = 'trial';
    end
  end

  states.uw_errors = sum(xor(tx_uw,rx_uw));

  % pre or post-amble has told us this is the start of the packet.  Confirm we 
  % have a valid frame by checking the UW after the modem frames containing
  % the UW have been received 
  if strcmp(states.sync_state,'trial')
    states.sync_counter++;
    if states.sync_counter == Nuwframes
      if states.uw_errors < states.bad_uw_errors;
        next_state = "synced";
        states.packet_count = 0;                          % number of packets in this burst
        states.modem_frame = Nuwframes;                   % which modem frame we are up to in packet
      else
        next_state = "search";
        % reset rxbuf to make sure we only ever do a postamble loop once through same samples
        states.rxbufst = states.Nrxbufhistory;
        states.rxbuf = zeros(1, states.Nrxbuf);
      end
    end
  end
  
  if strcmp(states.sync_state,'synced')
    states.modem_frame++;
    if (states.modem_frame >= states.Np) 
      states.modem_frame = 0;                           % start of new packet
      states.packet_count++;
      if (states.packetsperburst)
        if (states.packet_count >= states.packetsperburst)
          next_state = "search";                        % we've finished this burst
          % reset rxbuf to make sure we only ever do a postamble loop once through same samples
          states.rxbufst = states.Nrxbufhistory;
          states.rxbuf = zeros(1, states.Nrxbuf);
        end
      end
    end
  end
  
  states.last_sync_state = states.sync_state;
  states.sync_state = next_state;
endfunction

%-------------------------------------------------------
% fast sync voice state state_machine
%-------------------------------------------------------

function states = sync_state_machine_voice2(states, rx_uw)
  ofdm_load_const;
  next_state = states.sync_state;
  states.sync_start = states.sync_end = 0;

  if strcmp(states.sync_state,'search')

    if states.timing_valid
      states.frame_count = 0;
      states.sync_counter = 0;
      states.modem_frame = 0;
      states.sync_start = 1;
      next_state = 'trial';
    end
  end

  if strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial')

    states.frame_count++;

    % UW occurs at the start of a packet
    if states.modem_frame == 0
        states.uw_errors = sum(xor(tx_uw,rx_uw));

        if strcmp(states.sync_state,'trial')
          if states.uw_errors <= states.bad_uw_errors
            next_state = "synced";
          else
            next_state = "search";
          end
        end

        if strcmp(states.sync_state,'synced')
          if states.uw_errors > states.bad_uw_errors
            states.sync_counter++;
          else
            states.sync_counter = 0;
          end

          if states.sync_counter == 6
            next_state = "search";
          end
        end
      end

      % keep track of where we are up to in packet
      states.modem_frame++;
      if (states.modem_frame >= states.Np) states.modem_frame = 0; end
  end

  states.last_sync_state = states.sync_state;
  states.sync_state = next_state;
endfunction

