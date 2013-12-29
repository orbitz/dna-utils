let base = 5
let string_length = 1024

(* Uber lame *)
let rec pow b = function
  | 0 -> 1
  | n -> b * pow b (n - 1)

module Kmer = struct
  type t = { kmer     : int
	   ; pos      : int
	   ; max_pl   : int
	   ; kmer_rev : int array
	   }

  let int_of_nucleotide = function
    | 'a' | 'A' -> 0
    | 't' | 'T' -> 1
    | 'g' | 'G' -> 2
    | 'c' | 'C' -> 3
    | _         -> 4

  let nucleotide_of_int = function
    | 0 -> 'A'
    | 1 -> 'T'
    | 2 -> 'G'
    | 3 -> 'C'
    | 4 -> 'N'
    | _ -> failwith "impossible"

  let create kmer_len = { kmer     = 0
			; pos      = 0
			; max_pl   = pow base kmer_len
			; kmer_rev = Array.create kmer_len 0
			}

  let string_of_kmer kmer kmer_len =
    let s = String.create kmer_len in
    let rec sok' kmer = function
      | 0 -> s
      | n -> begin
	s.[n - 1] <- nucleotide_of_int (kmer mod base);
	sok' (kmer / base) (n - 1)
      end
    in
    sok' kmer kmer_len

  let add t c =
    let i        = int_of_nucleotide c in
    let l        = Array.length t.kmer_rev in
    let pos'     = (t.pos + 1) mod l in
    let rev_pos  = (t.pos + l) mod l in
    let kmer_rev = t.max_pl * t.kmer_rev.(rev_pos) in
    t.kmer_rev.(rev_pos) <- i;
    { t with
      kmer = t.kmer * base + i - kmer_rev
    ; pos  = pos'
    }

  let kmer t = t.kmer

end

type state = { data     : string
	     ; len      : int
	     ; kmers    : int array
	     ; kmer_len : int
	     }

let kmer_state kmer_len =
  let num_kmers = pow base kmer_len in
  { data     = String.create string_length
  ; len      = 0
  ; kmers    = Array.create num_kmers 0
  ; kmer_len = kmer_len
  }

let input_string ~pos ~len s =
  try
    input stdin s pos len
  with
    | End_of_file -> 0

let count_kmers ~data ~data_len ~kmers ~kmer_len =
  let rec init kmer kn n =
    if n < data_len then begin
      if kn < kmer_len - 1 then
	match data.[n] with
	  | '\n' | ' ' -> init kmer (kn + 1) (n + 1)
	  | '>'        -> skip_seq_header kmer n
	  | c          -> let kmer' = Kmer.add kmer c in
			  init kmer' (kn + 1) (n + 1)
      else
	proc kmer n
    end
    else
      n
  and proc kmer n =
    if n < data_len then begin
      match data.[n] with
	| '\n' | ' ' -> proc kmer (n + 1)
	| '>'        -> skip_seq_header kmer n
	| c          -> let kmer' = Kmer.add kmer c in
			kmers.(Kmer.kmer kmer') <- kmers.(Kmer.kmer kmer') + 1;
			proc kmer' (n + 1)
    end
    else
      n - kmer_len + 1
  and skip_seq_header kmer n =
    try
      match String.index_from data n '\n' with
      | n when n < data_len -> init (Kmer.create kmer_len) 0 (n + 1)
      | _                   -> n - kmer_len + 1
    with
      | Not_found -> n
  in
  init (Kmer.create kmer_len) 0 0

let rec loop = function
  | s when s.len < string_length -> begin
    let read =
      input_string
	~pos:s.len
	~len:(string_length - s.len)
	s.data
    in
    let s' = { s with len = s.len + read } in
    if read = 0 then begin
      ignore (count_kmers ~data:s'.data ~data_len:s'.len ~kmers:s'.kmers ~kmer_len:s'.kmer_len);
    end
    else
      loop s'
  end
  | s -> begin
    let processed =
      count_kmers
    	~data:s.data
    	~data_len:s.len
    	~kmers:s.kmers
    	~kmer_len:s.kmer_len
    in
    let len = s.len - processed in
    String.blit s.data processed s.data 0 len;
    loop { s with len = len }
  end


let main () =
  let kmer_len = int_of_string Sys.argv.(1) in
  (* This method can only handle kmer's upto 6 digits long *)
  assert (kmer_len > 0 && kmer_len <= 10);
  let s = kmer_state kmer_len in
  loop s;
  for i = 0 to Array.length s.kmers - 1 do
    if s.kmers.(i) > 0 then
      Printf.printf "%s\t%d\n" (Kmer.string_of_kmer i kmer_len) s.kmers.(i)
  done

let () = main ()
